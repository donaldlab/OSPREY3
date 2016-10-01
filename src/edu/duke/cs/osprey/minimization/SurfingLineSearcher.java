package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

public class SurfingLineSearcher implements LineSearcher {
	
	// NOTE: this class is based on the algorithm in the CCDMinimizer
	// this version is just easier to read and a little bit more optimized
	// TODO: port this to the GPU
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	public SurfingLineSearcher() {
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
	}
	
	@Override
	public void search(ObjectiveFunction f, DoubleMatrix1D x, int dof, DoubleMatrix1D mins, DoubleMatrix1D maxs) {
		
		double xdmin = mins.get(dof);
		double xdmax = maxs.get(dof);
		double xd = x.get(dof);
		
		double fxd = f.getValForDOF(dof, xd);
		Double fxdmin = null;
		Double fxdmax = null;
		
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
			fxdp = f.getValForDOF(dof, xdp);
		}
		double fxdm = Double.POSITIVE_INFINITY;
		if (xdm >= xdmin) {
			fxdm = f.getValForDOF(dof, xdm);
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

		// clamp xdstar to the range
		if (xdstar < xdmin) {
			xdstar = xdmin;
		}
		if (xdstar > xdmax) {
			xdstar = xdmax;
		}
		
		double fxdstar = f.getValForDOF(dof, xdstar);
		
		// did we go downhill?
		if (fxdstar < fxd) {
			
			// surf along f locally to try to find better minimum
			double xdsurfHere = xdstar;
			double fxdsurfHere = fxdstar;
			while (true) {
				
				// take a step twice as far as we did last time
				double xdsurfNext = xd + 2*(xdsurfHere - xd);
				
				// did we step off the min?
				if (xdsurfNext < xdmin) {
					
					// if the min is better, go there instead
					if (fxdmin == null) {
						fxdmin = f.getValForDOF(dof, xdmin);
					}
					if (fxdmin < fxdstar) {
						xdsurfHere = xdmin;
					}
					
					break;
				
				// did we step off the max?
				} else if (xdsurfNext > xdmax) {
					
					// if the max is better, go there instead
					if (fxdmax == null) {
						fxdmax = f.getValForDOF(dof, xdmax);
					}
					if (fxdmax < fxdstar) {
						xdsurfHere = xdmax;
					}
					
					break;
				}
				
				double fxdsurfNext = f.getValForDOF(dof, xdsurfNext);
				
				// did we improve the min enough to keep surfing?
				if (fxdsurfNext < fxdsurfHere - getTolerance(fxdsurfHere)) {
				
					// yeah, keep going
					xdsurfHere = xdsurfNext;
					fxdsurfHere = fxdsurfNext;
					
				} else {
					
					// nope, stop surfing
					break;
				}
			}
			
			// update the minimum estimate so far
			xdstar = xdsurfHere;
			fxdstar = fxdsurfHere;
			
		// did we go significantly uphill?
		} else if (fxdstar > fxd + Tolerance) {
			
			// try to surf back downhill
			double xdsurfHere = xdstar;
			double fxdsurfHere = fxdstar;
			while (true) {
				
				// cut the step in half
				double xdsurfNext = xd + (xdsurfHere - xd)/2;
				double fxdsurfNext = f.getValForDOF(dof, xdsurfNext);
				
				// did we improve the min enough to keep surfing?
				if (fxdsurfNext < fxdsurfHere - getTolerance(fxdsurfHere)) {
				
					// yeah, keep going
					xdsurfHere = xdsurfNext;
					fxdsurfHere = fxdsurfNext;
					
				} else {
					
					// nope, stop surfing
					break;
				}
			}
			
			// did the quadratic step help at all?
			if (fxdstar < fxd) {
				
				// yeah, keep it!
				
			} else {
				
				// nope, the original spot was lower
				xdstar = xd;
				fxdstar = fxd;
			}
			
			// did surfing help at all?
			if (fxdsurfHere < fxdstar) {
				
				// yeah, use the surf spot
				xdstar = xdsurfHere;
				fxdstar = fxdsurfHere;
			}
		}
		
		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum
		
		// NOTE: skipping this can make minimization a bit faster,
		// but skipping this causes a noticeable rise in final energies too
		// it's best to keep doing it I think
		
		xdm = xdstar - 1;
		if (xdm >= xdmin) {
			fxdm = f.getValForDOF(dof, xdm);
			if (fxdm < fxdstar) {
				xdstar = xdm;
				fxdstar = fxdm;
			}
		}
		
		xdp = xdstar + 1;
		if (xdp <= xdmax) {
			fxdp = f.getValForDOF(dof, xdp);
			if (fxdp < fxdstar) {
				xdstar = xdp;
				fxdstar = fxdp;
			}
		}
		
		lastStep = xdstar - xd;
		if (iteration == 0) {
			firstStep = lastStep;
		}

		// update the dofs and the conf
		x.set(dof, xdstar);
		f.setDOF(dof, xdstar);
		
		iteration++;
	}
	
	private double getTolerance(double f) {
		
		// use full tolerance, unless f is very small
		// then scale by the magnitude of f
		return Tolerance * Math.max(1, Math.abs(f));
	}
}
