package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.ojalgo.matrix.decomposition.LU;
import org.ojalgo.matrix.store.PrimitiveDenseStore;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.gpu.cuda.kernels.SubForcefieldsKernelCuda;
import edu.duke.cs.osprey.tools.Profiler;

public class QuasiNewtonMinimizer implements Minimizer {
	
	private static class Sample {
		
		public DoubleMatrix1D x;
		public double fx;
		
		public Sample(int n) {
			x = DoubleFactory1D.dense.make(n);
			fx = 0;
		}

		public void set(Sample other) {
			x.assign(other.x);
			fx = other.fx;
		}
	}
	
	private static final double MaxIterations = 30; // same as CCDMinimizer
	private static final double ConvergenceThreshold = 0.001; // same as CCDMinimizer
	
	private static final double Tolerance = 1e-6;
	private static final double InitialSurfaceDelta = 1; // for dihedral dofs, TODO: make configurable
	
	// TODO: play with step logic
	// looks like a bigger initial step size than 0.25 works better
	
	private ObjectiveFunction f;
	private ObjectiveFunction.DofBounds dofBounds;
	private List<ObjectiveFunction.OneDof> dofs;
	
	private SubForcefieldsKernelCuda kernel;
	
	public QuasiNewtonMinimizer(ObjectiveFunction f) {
		this(f, null);
	}
	
	public QuasiNewtonMinimizer(ObjectiveFunction f, SubForcefieldsKernelCuda kernel) {
		
		this.dofBounds = new ObjectiveFunction.DofBounds(f.getConstraints());
		
		// collect the the dofs
		this.f = f;
		this.dofs = new ArrayList<>();
		for (int d=0; d<f.getNumDOFs(); d++) {
			this.dofs.add(new ObjectiveFunction.OneDof(f, d));
		}
		
		this.kernel = kernel;
	}
	
	@Override
	public DoubleMatrix1D minimize() {
		
		int n = dofBounds.size();
		
		Profiler profiler = new Profiler();
		profiler.start("energy");
		
		// get the current objective function value
		// at the center of the bounds
		Sample here = new Sample(n);
		dofBounds.getCenter(here.x);
		here.fx = f.getValue(here.x);
		
		profiler.stop();
		
		DoubleMatrix1D dir = DoubleFactory1D.dense.make(n);
		Sample next = new Sample(n);
		
		double firstStep = InitialSurfaceDelta;
		double lastStep = InitialSurfaceDelta;
		
		for (int iter=0; iter<MaxIterations; iter++) {
			
			// update the surface delta, try to make it adaptive (based on historical steps if possible; else on step #)
			double surfaceDelta;
			if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
				surfaceDelta = InitialSurfaceDelta*Math.abs(lastStep/firstStep);
			} else {
				surfaceDelta = InitialSurfaceDelta/Math.pow(iter + 1, 3);
			}
			
			// TEMP
			//System.out.println(String.format("surface delta: %f   last: %f   prev: %f   ratio: %f", surfaceDelta, lastStep, prevStep, Math.abs(lastStep/prevStep)));
			
			profiler.resume("solve surface");
			Profiler stepProfiler = new Profiler();
			stepProfiler.start("solve surface");
			
			// solve for the local quadratic surface at x
			PrimitiveDenseStore surface = solveSurface(here, surfaceDelta);
			
			profiler.resume("minimize surface");
			stepProfiler.start("minimize surface");
			
			// analytically minimize the surface and get the direction relative to x
			DoubleMatrix1D xstar = minimizeSurfaceConstrained(surface);
			for (int d=0; d<n; d++) {
				dir.set(d, xstar.get(d) - here.x.get(d));
			}
			
			// TEMP
			//System.out.println(dump(dir));
			
			profiler.resume("line search");
			stepProfiler.start("line search");
			
			lineSearch(here, dir, next);
			
			profiler.stop();
			stepProfiler.stop();
			
			// TEMP
			//System.out.println(dump(dir));
			
			// update step metadata
			/* TODO: play with step choosing
			lastStep = 0;
			for (int d=0; d<n; d++) {
				lastStep = Math.max(lastStep, Math.abs(dir.get(d)));
			}
			*/
			lastStep = Math.sqrt(dir.zDotProduct(dir));
			if (iter == 0) {
				firstStep = lastStep;
			}
			
			// TEMP
			//System.out.println("effstep: " + lastStep);
			
			double improvement = here.fx - next.fx;
			
			// TEMP
			System.out.println(String.format("iter %3d   energy %12.6f   improvement %12.6f", iter, next.fx, improvement));
			//System.out.println(stepProfiler.makeReport(TimeUnit.MICROSECONDS));
			
			if (improvement > 0) {
				
				// keep the step
				here.set(next);
				
				if (improvement < ConvergenceThreshold) {
					break;
				}
				
			} else {
				break;
			}
		}
		
		// TEMP
		//System.out.println(profiler.makeReport(TimeUnit.MILLISECONDS));
		
		// update the protein
		f.setDOFs(here.x);
		
		return here.x;
	}
	
	private String dump(DoubleMatrix1D x) {
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<x.size(); i++) {
			if (i > 0) {
				buf.append(" ");
			}
			buf.append(String.format("%8.3f", x.get(i)));
		}
		return buf.toString();
	}
	
	private String dump(PrimitiveDenseStore x) {
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<x.count(); i++) {
			if (i > 0) {
				buf.append(" ");
			}
			buf.append(String.format(" %8.3f", x.get(i)));
		}
		return buf.toString();
	}
	
	private PrimitiveDenseStore solveSurface(Sample here, double delta) {
		
		int n = dofBounds.size();
		int m = 2*n + 1;
		PrimitiveDenseStore A = PrimitiveDenseStore.FACTORY.makeZero(m, m);
		PrimitiveDenseStore b = PrimitiveDenseStore.FACTORY.makeZero(m, 1);
		
		int r = 0;
		DoubleMatrix1D y = DoubleFactory1D.dense.make(n);
		
		// sample the center point
		setSurfaceRow(A, b, r++, here.x, here.fx);
		
		// do we have a gpu kernel?
		if (kernel != null) {
			
			// yup, compute the dof energies in parallel
			f.setDOFs(here.x);
			kernel.calcEnergies(here.x, delta);
			
			// copy the samples into A,b
			for (int d=0; d<n; d++) {
				
				y.assign(here.x);
				double xd = here.x.get(d);
				
				y.set(d, xd - delta);
				setSurfaceRow(A, b, r++, y, here.fx + kernel.getFxdmOffset(d));
				
				y.set(d, xd + delta);
				setSurfaceRow(A, b, r++, y, here.fx + kernel.getFxdpOffset(d));
				
				// TODO: make a unit test to check for accuracy
				//System.out.println(String.format("\td %2d   m %12.6f   p %12.6f", d, kernel.getFxdmOffset(d), kernel.getFxdpOffset(d)));
			}
			
		} else {
			
			// nope, compute the dof energies in serial
			for (int d=0; d<n; d++) {
				
				// reset the protein to x
				f.setDOFs(here.x);
				y.assign(here.x);
				
				double xd = here.x.get(d);
				ObjectiveFunction.OneDof dof = dofs.get(d);
				
				double fxd = dof.getValue(xd);
				
				// sample the minus delta point
				double xdm = xd - delta;
				y.set(d, xdm);
				double fxdm = dof.getValue(xdm);
				setSurfaceRow(A, b, r++, y, here.fx - fxd + fxdm);
				
				// sample the plus delta point
				double xdp = xd + delta;
				y.set(d, xdp);
				double fxdp = dof.getValue(xdp);
				setSurfaceRow(A, b, r++, y, here.fx - fxd + fxdp);
				
				// TEMP
				//System.out.println(String.format("\td %2d   m %12.6f   p %12.6f", d, fxdm - fxd, fxdp - fxd));
			}
		}
		
		assert (r == m);
		
		// solve the linear system Ax=b
		LU<Double> solver = LU.PRIMITIVE.make();
		solver.decompose(A);
		assert (solver.isSolvable());
		PrimitiveDenseStore surface = (PrimitiveDenseStore)solver.solve(b);
		
		/* DEBUG: check the surface for accuracy
		checkSurface(0, here.x, y, 0, surface);
		for (int d=0; d<n; d++) {
			checkSurface(d, here.x, y, -delta, surface);
			checkSurface(d, here.x, y, delta, surface);
		}
		*/
		
		return surface;
	}
	
	private void checkSurface(int d, DoubleMatrix1D x, DoubleMatrix1D y, double delta, PrimitiveDenseStore surface) {
		
		// update y
		y.assign(x);
		y.set(d, y.get(d) + delta);
		
		// get the real energy
		double fy = f.getValue(y);
		
		// evaluate the quadratic surface
		int r = 0;
		double qy = surface.get(r++, 0);
		for (int d2=0; d2<dofBounds.size(); d2++) {
			double yd = y.get(d2);
			qy += surface.get(r++, 0)*yd;
			qy += surface.get(r++, 0)*yd*yd;
		}
		
		System.out.println(String.format("\td: %2d   delta: %7.4f   fx: %12.6f   fqx: %12.6f   err: %.12f", d, delta, fy, qy, Math.abs(fy - qy)));
	}

	private void setSurfaceRow(PrimitiveDenseStore A, PrimitiveDenseStore b, int r, DoubleMatrix1D x, double fx) {
	
		int n = dofBounds.size();
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
	
	// sadly, this doesn't work as well as using the constraints
	// it's harder to pick a step size when the goal could be outside the bounds
	private DoubleMatrix1D minimizeSurfaceUnconstrained(PrimitiveDenseStore surface) {
		
		// minimize the quadratic surface (but ignore the bounds)
		// convert it to a system of linear equations to solve for the vertex
		// luckily, all the equations are independent, so we can use regular algebra
		// otherwise, we could use linear algebra (like LU decomposition)
		
		// NOTE: I tried including cross-dimensional terms in the quadratic surface,
		// but this resulted in poor minimization quality (ie, higher energies)
		
		int n = dofBounds.size();
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
		
		int n = dofBounds.size();
		
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
				return dofBounds.getBounds();
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
				return InitialSurfaceDelta;
			}

			@Override
			public boolean isDOFAngle(int dof) {
				return true;
			}
		};
		
		// minimize it!
		return new CCDMinimizer(fq, false).minimize();
	}
	
	public void lineSearch(Sample here, DoubleMatrix1D dir, Sample out) {
		
		int n = here.x.size();
		
		// see where a full step would take us
		double step = 1;
		takeStep(out, here, dir, step);
		out.fx = f.getValue(out.x);
		
		// surf!
		double stepNext = 0;
		Sample next = new Sample(n);
		next.set(out);
		while (true) {
			
			// cut the step in half
			stepNext = step/2;
			takeStep(next, here, dir, stepNext);
			next.fx = f.getValue(next.x);
			
			// did we improve enough to keep surfing?
			if (next.fx < out.fx - getTolerance(out.fx)) {
			
				// yeah, keep going
				step = stepNext;
				out.set(next);
				
			} else {
				
				// nope, stop surfing
				break;
			}
		}
		
		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum
		
		takeStep(next, here, dir, stepNext, -1);
		if (dofBounds.isInBounds(next.x)) {
			next.fx = f.getValue(next.x);
			if (next.fx < out.fx) {
				out.set(next);
			}
		}
		
		takeStep(next, here, dir, stepNext, 1);
		if (dofBounds.isInBounds(next.x)) {
			next.fx = f.getValue(next.x);
			if (next.fx < out.fx) {
				out.set(next);
			}
		}
		
		// update dir
		for (int d=0; d<n; d++) {
			dir.set(d, out.x.get(d) - here.x.get(d));
		}
	}

	private void takeStep(Sample out, Sample here, DoubleMatrix1D dir, double step) {
		takeStep(out, here, dir, step, 0);
	}
	
	private void takeStep(Sample out, Sample here, DoubleMatrix1D dir, double step, double offset) {
		
		int n = here.x.size();
		Double length = null;
		
		for (int d=0; d<n; d++) {
			double xd = here.x.get(d);
			
			xd += dir.get(d)*step;
			
			if (offset != 0) {
				if (length == null) {
					length = Math.sqrt(dir.zDotProduct(dir));
				}
				xd += dir.get(d)*offset/length;
			}
			
			out.x.set(d, xd);
		}
	}
	
	private double getTolerance(double f) {
		
		// use full tolerance, unless f is very small
		// then scale by the magnitude of f
		return Tolerance * Math.max(1, Math.abs(f));
	}
}
