package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.ojalgo.matrix.decomposition.LU;
import org.ojalgo.matrix.store.PrimitiveDenseStore;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.tools.Profiler;

public class QuasiNewtonMinimizer implements Minimizer.NeedsCleanup {
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private ObjectiveFunction f;
	private List<ObjectiveFunction.OneDof> dofs;
	
	private boolean constrainSurface;
	
	public QuasiNewtonMinimizer(ObjectiveFunction f) {
		
		// using the constraints on the quadratic surface minimization seems to give better accuracy
		// so it should probably be on by default
		// turning off constraints gives significantly faster performance, but accuracy suffers a bit
		this(f, true);
	}

	public QuasiNewtonMinimizer(ObjectiveFunction f, boolean constrainSurface) {
		
		this.f = f;
		this.constrainSurface = constrainSurface;
		
		// collect the the dofs
		dofs = new ArrayList<>();
		for (int d=0; d<f.getNumDOFs(); d++) {
			ObjectiveFunction.OneDof fd = new ObjectiveFunction.OneDof(f, d);
			dofs.add(fd);
		}
	}
	
	@Override
	public DoubleMatrix1D minimize() {
		
		Profiler profiler = new Profiler();
		
		// init x to the center of the bounds
		int n = dofs.size();
		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);
		for (int d=0; d<n; d++) {
			ObjectiveFunction.OneDof dof = dofs.get(d);
			x.set(d, (dof.getXMin() + dof.getXMax())/2);
		}
		
		DoubleMatrix1D dir = DoubleFactory1D.dense.make(n);
		
		profiler.resume("energy");
		
		// get the current objective function value
		double curf = f.getValue(x);
		
		profiler.stop();
		
		double firstStep = 1;
		double lastStep = 1;
		
		for (int iter=0; iter<MaxIterations; iter++) {
			
			// update the surface delta, try to make it adaptive (based on historical steps if possible; else on step #)
			double surfaceDelta;
			if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
				surfaceDelta = InitialStepSize*Math.abs(lastStep/firstStep);
			} else {
				surfaceDelta = InitialStepSize/Math.pow(iter + 1, 3);
			}
			
			// TEMP
			//System.out.println(String.format("surface delta: %f", surfaceDelta));
			
			profiler.resume("solve surface");
			
			// solve for the local quadratic surface at x
			PrimitiveDenseStore surface = solveSurface(x, surfaceDelta, curf);
			
			if (constrainSurface) {
				
				profiler.resume("minimize surface");
				
				// analytically minimize the surface and get the direction relative to x
				DoubleMatrix1D xstar = minimizeSurface(surface, true);
				for (int d=0; d<n; d++) {
					dir.set(d, xstar.get(d) - x.get(d));
				}
				
				profiler.resume("line search");
				
				lineSearch(x, dir);
				
				profiler.stop();
				
			} else {
				
				profiler.resume("minimize surface");
				
				// analytically minimize the surface and get the direction relative to x
				DoubleMatrix1D xstar = minimizeSurface(surface, false);
				for (int d=0; d<n; d++) {
					dir.set(d, xstar.get(d) - x.get(d));
				}
				
				// TEMP
				//System.out.println("direction: " + dump(dir));
				
				// damp the direction to avoid weighting steep dimensions too heavily
				// (helps to avoid being too "scared" of clashes)
					
				// magic number, chosen empirically based on limited test cases
				// seems to work well in practice, but could be sub-optimal in some cases
				final double DampRate = 1.1;
				
				for (int d=0; d<n; d++) {
					double xd = dir.get(d);
					
					if (xd > 0) {
						xd = Math.log(xd + 1)*DampRate;
					} else if (xd < 0) {
						xd = -Math.log(1 - xd)*DampRate;
					}
					
					// NOTE: this damping function seems to work well too,
					// but not as well as the log
					//xd = Math.atan(xd);
					
					dir.set(d, xd);
				}
				
				// TEMP
				//System.out.println("damped dr: " + dump(dir));
				
				// hard clamp the dir so we stay within the bounds
				clampDir(x, dir);
				
				/* NOTE: don't try to preserve the step direction with scaling.
					this causes the steps to heavily favor steeper dimensions.
					attempting to account for that by ignoring dimensions that hit a bound
					doesn't perform as well as this simpler damp-and-clamp method
				*/
				
				profiler.stop();
			}
			
			profiler.resume("step");
			
			// update step metadata
			lastStep = Math.sqrt(dir.zDotProduct(dir));
			if (iter == 0) {
				firstStep = lastStep;
			}
			
			// TEMP
			//System.out.println("effstep: " + lastStep);
			
			// take the step
			for (int d=0; d<n; d++) {
				x.set(d, x.get(d) + dir.get(d));
			}
			
			profiler.resume("energy");
			
			double nextf = f.getValue(x);
			
			profiler.stop();
			
			// TEMP
			System.out.println(String.format("iter %3d   energy %12.6f   improvement %12.6f", iter, nextf, curf - nextf));
			
			// did we improve enough to keep going?
			if (curf - nextf < ConvergenceThreshold) {
				
				// nope, we're done
				break;
				
			} else {
				
				// yeah, keep going
				curf = nextf;
			}
		}
		
		System.out.println(profiler.makeReport(TimeUnit.MILLISECONDS));
		
		return x;
	}
	
	private String dump(DoubleMatrix1D x) {
		StringBuilder buf = new StringBuilder();
		int n = x.size();
		for (int d=0; d<n; d++) {
			if (d > 0) {
				buf.append(" ");
			}
			buf.append(String.format("%8.3f", x.get(d)));
		}
		return buf.toString();
	}

	private PrimitiveDenseStore solveSurface(DoubleMatrix1D x, double delta, double fx) {
		
		int n = dofs.size();
		int m = 2*n + 1;
		PrimitiveDenseStore A = PrimitiveDenseStore.FACTORY.makeZero(m, m);
		PrimitiveDenseStore b = PrimitiveDenseStore.FACTORY.makeZero(m, 1);
		
		int r = 0;
		
		// sample the center point
		setSurfaceRow(A, b, r++, x, fx);
		
		DoubleMatrix1D y = x.copy();
		for (int d=0; d<n; d++) {
			
			y.assign(x);
			
			double xd = x.get(d);
			ObjectiveFunction.OneDof dof = dofs.get(d);
			
			double fxd = dof.getValue(xd);
			
			// sample the minus delta point
			double xdm = xd - delta;
			y.set(d, xdm);
			double fxdm = dof.getValue(xdm);
			//double fxdm = f.getValue(y);
			//setSurfaceRow(numActiveDofs, A, b, r++, y, fxdm);
			setSurfaceRow(A, b, r++, y, fx - fxd + fxdm);
			
			// sample the plus delta point
			double xdp = xd + delta;
			y.set(d, xdp);
			double fxdp = dof.getValue(xdp);
			//double fxdp = f.getValue(y);
			//setSurfaceRow(numActiveDofs, A, b, r++, y, fxdp);
			setSurfaceRow(A, b, r++, y, fx - fxd + fxdp);
		}
		
		assert (r == m);
		
		// solve the linear system Ax=b
		LU<Double> solver = LU.PRIMITIVE.make(A);
		solver.decompose(A);
		assert (solver.isSolvable());
		return (PrimitiveDenseStore)solver.solve(b);
	}

	private void setSurfaceRow(PrimitiveDenseStore A, PrimitiveDenseStore b, int r, DoubleMatrix1D x, double fx) {
	
		int n = dofs.size();
		int c = 0;
		
		// constant term
		A.set(r, c++, 1);
		
		for (int d=0; d<n; d++) {
			double xd = x.get(d);
			
			// linear terms
			A.set(r, c++, xd);
			
			// quadratic terms
			A.set(r, c++, xd*xd);
		}
		
		assert (c == 2*n + 1);
		
		b.set(r, fx);
	}
	
	private DoubleMatrix1D minimizeSurface(PrimitiveDenseStore surface, boolean useConstraints) {
		if (useConstraints) {
			return minimizeSurfaceConstrained(surface);
		} else {
			return minimizeSurfaceUnconstrained(surface);
		}
	}
	
	private DoubleMatrix1D minimizeSurfaceUnconstrained(PrimitiveDenseStore surface) {
		
		// minimize the quadratic surface (but ignore the bounds)
		// convert it to a system of linear equations to solve for the vertex
		// luckily, all the equations are independent, so we can use regular algebra
		// otherwise, we could use linear algebra (like LU decomposition)
		
		// NOTE: I tried including cross-dimensional terms in the quadratic surface,
		// but this resulted in poor minimization quality (ie, higher energies)
		
		int n = dofs.size();
		DoubleMatrix1D x = DoubleFactory1D.dense.make(n);
		
		// ignore the first row, it has the constant value which we don't need anyway
		int r = 1;
		
		for (int d=0; d<n; d++) {
			
			// solve for the minimal xd
			// ax^2 + bx = 0
			double b = surface.get(r++, 0);
			double a = surface.get(r++, 0);
			double xdstar = -b/2/a;
			
			x.set(d, xdstar);
		}
		
		return x;
	}
	
	private DoubleMatrix1D minimizeSurfaceConstrained(PrimitiveDenseStore surface) {
		
		// minimize the constrained quadratic surface
		// ie don't necessarily choose the vertex
		// TODO: use a QP solver. it'll probably be faster than CCD
		
		int n = dofs.size();
		
		// make an objective function for the quadratic surface
		ObjectiveFunction fq = new ObjectiveFunction() {
			
			private static final long serialVersionUID = 4038936433792726958L;
			
			private DoubleMatrix1D x;
			private double constant;
			private double[] linears;
			private double[] quadratics;
			
			{
				// allocate space
				x = DoubleFactory1D.dense.make(n);
				linears = new double[n];
				quadratics = new double[n];
				
				// unpack the quadratic surface coefficients
				int r = 0;
				constant = surface.get(r++, 0);
				for (int d=0; d<n; d++) {
					linears[d] = surface.get(r++, 0);
					quadratics[d] = surface.get(r++, 0);
				}
				assert (r == 2*n + 1);
			}
			
			@Override
			public int getNumDOFs() {
				return n;
			}

			@Override
			public DoubleMatrix1D[] getConstraints() {
				return f.getConstraints();
			}

			@Override
			public void setDOFs(DoubleMatrix1D x) {
				this.x.assign(x);
			}

			@Override
			public void setDOF(int d, double xd) {
				x.set(d, xd);
			}

			@Override
			public double getValue(DoubleMatrix1D x) {
				setDOFs(x);
				
				double val = constant;
				for (int d=0; d<n; d++) {
					double xd = this.x.get(d);
					val += linears[d]*xd;
					val += quadratics[d]*xd*xd;
				}
				
				return val;
			}

			@Override
			public double getValForDOF(int d, double xd) {
				setDOF(d, xd);
				
				return constant + linears[d]*xd + quadratics[d]*xd*xd;
			}

			@Override
			public double getInitStepSize(int dof) {
				return InitialStepSize;
			}

			@Override
			public boolean isDOFAngle(int dof) {
				return true;
			}
		};
		
		// minimize it!
		return new CCDMinimizer(fq, false).minimize();
	}
	
	private void clampDir(DoubleMatrix1D x, DoubleMatrix1D dir) {
		
		int n = x.size();
		
		for (int d=0; d<n; d++) {
			
			double xd = x.get(d);
			double yd = xd + dir.get(d);
			
			ObjectiveFunction.OneDof dof = dofs.get(d);
			if (yd < dof.getXMin()) {
				dir.set(d, dof.getXMin() - xd);
			} else if (yd > dof.getXMax()) {
				dir.set(d, dof.getXMax() - xd);
			}
		}
	}

	private boolean isInBounds(DoubleMatrix1D x) {
		int n = x.size();
		for (int d=0; d<n; d++) {
			double xd = x.get(d);
			
			ObjectiveFunction.OneDof dof = dofs.get(d);
			if (xd < dof.getXMin() || xd > dof.getXMax()) {
				return false;
			}
		}
		
		return true;
	}

	public void lineSearch(DoubleMatrix1D x, DoubleMatrix1D dir) {
		
		int n = x.size();
		
		// see where a full step would take us
		double step = 1;
		DoubleMatrix1D y = DoubleFactory1D.dense.make(n);
		updateY(y, x, dir, step);
		double fy = f.getValue(y);
		
		DoubleMatrix1D yNext = DoubleFactory1D.dense.make(n); 
		double fyNext = 0;
		double stepNext = 0;
		
		// surf!
		while (true) {
			
			// cut the step in half
			stepNext = step/2;
			
			// check f
			updateY(yNext, x, dir, stepNext);
			fyNext = f.getValue(yNext);
			
			// did we improve enough to keep surfing?
			if (fyNext < fy - getTolerance(fy)) {
			
				// yeah, keep going
				step = stepNext;
				y.assign(yNext);
				fy = fyNext;
				
			} else {
				
				// nope, stop surfing
				break;
			}
		}
		
		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum
		
		updateY(yNext, x, dir, stepNext, -1);
		if (isInBounds(yNext)) {
			fyNext = f.getValue(yNext);
			if (fyNext < fy) {
				y.assign(yNext);
				fy = fyNext;
			}
		}
		
		updateY(yNext, x, dir, stepNext, 1);
		if (isInBounds(yNext)) {
			fyNext = f.getValue(yNext);
			if (fyNext < fy) {
				y.assign(yNext);
				fy = fyNext;
			}
		}
		
		// update dir
		for (int d=0; d<n; d++) {
			dir.set(d, y.get(d) - x.get(d));
		}
	}

	private void updateY(DoubleMatrix1D y, DoubleMatrix1D x, DoubleMatrix1D dir, double step) {
		updateY(y, x, dir, step, 0);
	}
	
	private void updateY(DoubleMatrix1D y, DoubleMatrix1D x, DoubleMatrix1D dir, double step, double offset) {
		
		int n = x.size();
		Double length = null;
		
		for (int d=0; d<n; d++) {
			double xd = x.get(d);
			
			xd += dir.get(d)*step;
			
			if (offset != 0) {
				if (length == null) {
					length = Math.sqrt(dir.zDotProduct(dir));
				}
				xd += dir.get(d)*offset/length;
			}
			
			y.set(d, xd);
		}
	}
	
	private double getTolerance(double f) {
		
		// use full tolerance, unless f is very small
		// then scale by the magnitude of f
		return Tolerance * Math.max(1, Math.abs(f));
	}

	@Override
	public void cleanup() {
		// TODO
	}
}
