package edu.duke.cs.osprey.minimization;

import cern.colt.matrix.DoubleMatrix1D;

public class GpuStyleSurfingLineSearcher implements LineSearcher {
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	private ObjectiveFunction f;
	private int dof;
	private double xdmin;
	private double xdmax;
	private double xd;
	private double fxd;
	private double fxdmin;
	private double fxdmax;
	private double step;
	private double xdm;
	private double xdp;
	private double fxdm;
	private double fxdp;
	private double xdstar;
	private double fxdstar;
	
	private double stepScale;
	private boolean isSurfing;
	private int surfStep;
	private double xdsurfHere;
	private double fxdsurfHere;
	
	private boolean doEnergy;
	private double energyxd;
	private double energyfxd;
	
	public GpuStyleSurfingLineSearcher() {
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
	}
	
	@Override
	public void search(ObjectiveFunction f, DoubleMatrix1D x, int dof, DoubleMatrix1D mins, DoubleMatrix1D maxs) {
		
		this.f = f;
		this.dof = dof;
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// pretend like each function call is a kernel launch
		// so they can only communicate via shared memory
		
		// init memory
		xdmin = mins.get(dof);
		xdmax = maxs.get(dof);
		xd = x.get(dof);
		
		// init first energy calc
		doEnergy = true;
		energyxd = xd;
		
		// INPUTS: xdmin, xdmax, xd, doEnergy, energyxd
		
		energy();
		search1();
		energy();
		search2();
		energy();
		search3();
		energy();
		search4();
		energy();
		search5();
		energy();
		search6();
		for (int i=0; i<10; i++) {
			searchSurf();
			energy();
		}
		search7();
		energy();
		search8();
		energy();
		search9();
		
		// OUTPUTS: step, xdstar
		
		// update the step
		lastStep = step;
		if (iteration == 0) {
			firstStep = lastStep;
		}

		// update the dofs and the conf
		x.set(dof, xdstar);
		f.setDOF(dof, xdstar);
		
		iteration++;
	}
	
	private void energy() {
		if (doEnergy) {
			energyfxd = f.getValForDOF(dof, energyxd);
		}
	}
	
	private void search1() {
		
		// read energy result
		fxd = energyfxd;
		
		// get the min edge
		doEnergy = true;
		energyxd = xdmin;
	}
	
	private void search2() {
		
		// read energy result
		fxdmin = energyfxd;
		
		// get the max edge
		doEnergy = true;
		energyxd = xdmax;
	}
	
	private void search3() {
		
		// read energy result
		fxdmax = energyfxd;
		
		// get the minus neighbor
		xdm = xd - step;
		
		// init next energy calc
		doEnergy = xdm >= xdmin;
		energyxd = xdm;
	}
	
	private void search4() {
		
		// read energy result
		fxdm = Double.POSITIVE_INFINITY;
		if (doEnergy) {
			fxdm = energyfxd;
		}
		
		// get the plus neighbor
		xdp = xd + step;
		
		// init next energy calc
		doEnergy = xdp <= xdmax;
		energyxd = xdp;
	}
	
	private void search5() {
		
		// read energy result
		fxdp = Double.POSITIVE_INFINITY;
		if (doEnergy) {
			fxdp = energyfxd;
		}
		
		// fit a quadratic to the objective function, locally:
		// q(x) = fx + a*(x - xd)^2 + b*(x - xd)
		// a*step^2 + b*step = fxp - fx
		// a*step^2 - b*step = fxm - fx
		
		// solve for a to determine the shape
		double a = (fxdp + fxdm - 2*fxd)/(2*step*step);
		xdstar = 0;
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
		
		// init next energy calc
		doEnergy = true;
		energyxd = xdstar;
	}
	
	private void search6() {
		
		// read energy result
		fxdstar = energyfxd;
		
		// if we went downhill, use growing steps
		if (fxdstar < fxd) {
			stepScale = 2;
		} else {
			stepScale = 0.5;
		}
		
		// surf along f locally to try to find better minimum
		isSurfing = true;
		surfStep = 0;
		xdsurfHere = xdstar;
		fxdsurfHere = fxdstar;
	}
	
	private void searchSurf() {
		
		doEnergy = false;
			
		if (isSurfing && surfStep > 0) {
				
			double xdsurfNext = energyxd;
			double fxdsurfNext = energyfxd;
			
			// did we improve the min enough to keep surfing?
			if (fxdsurfNext < fxdsurfHere - getTolerance(fxdsurfHere)) {
			
				// yeah, keep going
				xdsurfHere = xdsurfNext;
				fxdsurfHere = fxdsurfNext;
				
			} else {
				
				// nope, stop surfing
				isSurfing = false;
			}
		}
			
		if (isSurfing) {
		
			// surf a step
			double xdsurfNext = xd + stepScale*(xdsurfHere - xd);
			
			// did we step off the min?
			if (xdsurfNext < xdmin) {
				
				// if the min is better, go there instead
				if (fxdmin < fxdstar) {
					xdsurfHere = xdmin;
				}
				
				isSurfing = false;
			
			// did we step off the max?
			} else if (xdsurfNext > xdmax) {
				
				// if the max is better, go there instead
				if (fxdmax < fxdstar) {
					xdsurfHere = xdmax;
				}
				
				isSurfing = false;
				
			} else {
			
				// init energy calc for new surf spot
				surfStep++;
				doEnergy = true;
				energyxd = xdsurfNext;
			}
		}
	}
	
	private void search7() {
		
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
		
		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum
		
		// NOTE: skipping this can make minimization a bit faster,
		// but skipping this causes a noticeable rise in final energies too
		// it's best to keep doing it I think
		
		// init next energy calc
		double xdm = xdstar - 1;
		doEnergy = xdm >= xdmin;
		energyxd = xdm;
	}
	
	private void search8() {
		
		// read energy result
		if (doEnergy) {
			if (energyfxd < fxdstar) {
				xdstar = energyxd;
				fxdstar = energyfxd;
			}
		}
		
		// init next energy calc
		double xdp = xdstar + 1;
		doEnergy = xdp <= xdmax;
		energyxd = xdp;
	}
	
	private void search9() {
		
		// read energy result
		if (doEnergy) {
			if (energyfxd < fxdstar) {
				xdstar = energyxd;
				fxdstar = energyfxd;
			}
		}
		
		// update the step to the one we actually took
		step = xdstar - xd;
	}
	
	private double getTolerance(double f) {
		
		// use full tolerance, unless f is very small
		// then scale by the magnitude of f
		return Tolerance * Math.max(1, Math.abs(f));
	}
}
