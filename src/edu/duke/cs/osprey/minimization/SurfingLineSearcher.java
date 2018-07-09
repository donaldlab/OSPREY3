/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.minimization;

public class SurfingLineSearcher implements LineSearcher {
	
	// NOTE: this class is based on the algorithm in the CCDMinimizer
	// this version is just easier to read and a little bit more optimized
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private ObjectiveFunction.OneDof f;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	public SurfingLineSearcher() {
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
		Double fxdmin = null;
		Double fxdmax = null;
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		double step;
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// make sure the step isn't so big that the quadratic approximation is worthless
		while (xd - step < xdmin && xd + step > xdmax) {
			step /= 2;
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
		
		// solve for the shape of the parabola
		double shape = fxdp + fxdm - 2*fxd;
		final double ShapeEpsilon = 1e-12;
		double xdstar = 0;
		if ((shape < -ShapeEpsilon) || Double.isNaN(shape) || Double.isInfinite(shape)) {
			
			// negative shape means quadratic is concave down
			// infinite or nan a means we're hitting a constraint or impossible conformation
			// so just minimize over the endpoints of the interval
			if (fxdm < fxdp) {
				xdstar = xdm;
			} else {
				xdstar = xdp;
			}
			
		} else if (shape <= ShapeEpsilon) {
			
			// shape near zero means it's basically flat here
			// so don't step anywhere
			xdstar = xd;
			
		} else {
			
			// positive shape means quadratic is concave up
			// step to the optimum
			
			/* this isn't terribly numerically stable, don't use it
			// solve for a and b params
			double a = shape/(2*step*step);
			double b = (fxdp - fxd)/step - a*step;
			// then minimize the quadratic to get the minimum x:
			// 2*a*(x - xd) + b = 0
			double deltax = -b/2/a;
			*/
			
			// this is much more numerically stable
			xdstar = xd + (fxdm - fxdp)*step/2/shape;
		}
		
		// clamp xdstar to the range
		if (xdstar < xdmin) {
			xdstar = xdmin;
		}
		if (xdstar > xdmax) {
			xdstar = xdmax;
		}
		
		double fxdstar = f.getValue(xdstar);
		
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
						fxdmin = f.getValue(xdmin);
					}
					if (fxdmin < fxdsurfHere) {
						xdsurfHere = xdmin;
						fxdsurfHere = fxdmin;
					}
					
					break;
				
				// did we step off the max?
				} else if (xdsurfNext > xdmax) {
					
					// if the max is better, go there instead
					if (fxdmax == null) {
						fxdmax = f.getValue(xdmax);
					}
					if (fxdmax < fxdsurfHere) {
						xdsurfHere = xdmax;
						fxdsurfHere = fxdmax;
					}
					
					break;
				}
				
				double fxdsurfNext = f.getValue(xdsurfNext);
				
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
				double fxdsurfNext = f.getValue(xdsurfNext);
				
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
		
		// update step before wall jumping
		lastStep = xdstar - xd;
		if (iteration == 0) {
			firstStep = lastStep;
		}
		
		// try to jump over walls arbitrarily
		// look in a 1-degree step for a better minimum
		
		// NOTE: skipping this can make minimization a bit faster,
		// but skipping this causes a noticeable rise in final energies too
		// it's best to keep doing it I think
		
		xdm = xdstar - 1;
		xdp = xdstar + 1;
		
		if (xdm >= xdmin) {
			fxdm = f.getValue(xdm);
			if (fxdm < fxdstar) {
				xdstar = xdm;
				fxdstar = fxdm;
			}
		}
		
		if (xdp <= xdmax) {
			fxdp = f.getValue(xdp);
			if (fxdp < fxdstar) {
				xdstar = xdp;
				fxdstar = fxdp;
			}
		}
		
		iteration++;
		
		f.setX(xdstar);
		return xdstar;
	}
	
	private double getTolerance(double f) {
		
		// scale abs(f) by tolerance, unless f is very small
		return Tolerance * Math.max(1, Math.abs(f));
	}
}
