package edu.duke.cs.osprey.minimization;

public class QuadraticLineSearcher implements LineSearcher {
	
	// NOTE: this class is based on the algorithm in the CCDMinimizer
	// this version is just easier to read and a little bit more optimized
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private ObjectiveFunction.OneDof f;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	public QuadraticLineSearcher() {
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
	}
	
	@Override
	public void init(ObjectiveFunction.OneDof f) {
		this.f = f;
	}
	
	@Override
	public double search(double xd) {
		
		double xdmin = f.getXMin();
		double xdmax = f.getXMax();
		
		double fxd = f.getValue(xd);
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		double step;
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// get the positive (p) and negative (n) neighbors for our current pos
		double xdp = xd + step;
		double xdm = xd - step;
		double fxdp = Double.POSITIVE_INFINITY;
		if (xdp <= xdmax) {
			fxdp = f.getValue(xdp);
		}
		double fxdm = Double.POSITIVE_INFINITY;
		if (xdm >= xdmin) {
			fxdm = f.getValue(xdm);
		}
		
		// fit a quadratic to the objective function, locally:
		// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
		// a*step^2 + b*step = fxp - fx
		// a*step^2 - b*step = fxm - fx
		
		// solve for a to determine the shape
		double a = (fxdp + fxdm - 2*fxd)/(2*step*step);
		double xdstar = 0;
		if ((a <= 0) || Double.isNaN(a) || Double.isInfinite(a)) {
			
			// negative a means quadratic is concave down, I think
			// infinite or nan a means we're hitting a constraint or impossible conformation
			// so just minimize over the endpoints of the interval
			if (fxdm < fxdp) {
				xdstar = xdm;
			} else {
				xdstar = xdp;
			}
			
		} else {
			
			// positive a means quadratic is concave up, I think
			// solve for the b param
			double b = (fxdp - fxd)/step - a*step;
			
			// then minimize the quadratic to get the minimum x:
			// 2*a*(x - xd) + b = 0
			
			xdstar = xd - b/2/a;
		}
		
		/* TEMP
		// clamp xdstar to the range
		if (xdstar < xdmin) {
			xdstar = xdmin;
		}
		if (xdstar > xdmax) {
			xdstar = xdmax;
		}
		*/
		
		lastStep = xdstar - xd;
		if (iteration == 0) {
			firstStep = lastStep;
		}
		
		iteration++;
		
		f.setX(xdstar);
		return xdstar;
	}
}
