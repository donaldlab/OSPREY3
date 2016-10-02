package edu.duke.cs.osprey.minimization;

import static edu.duke.cs.osprey.tools.VectorAlgebra.*;

import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Residue;

public class GpuStyleSurfingLineSearcher implements LineSearcher {
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs, TODO: make configurable
	
	private EnergyFunction efunc;
	private FreeDihedral dof;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
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
	
	private double[] coords;
	private int[] dihedralAtomIndices;
	private int[] rotatedAtomIndices;
	
	private boolean doEnergy;
	private double energyxd;
	private double energyfxd;
	
	public GpuStyleSurfingLineSearcher() {
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
	}
	
	@Override
	public void search(ObjectiveFunction f, DoubleMatrix1D x, int d, DoubleMatrix1D mins, DoubleMatrix1D maxs) {
		
		// HACKHACK: break the interfaces to get the pieces we need
		// the runtime will complain if we can't do these casts, so we're still safe, but the compiler has no idea what's going on
		MoleculeModifierAndScorer mof = (MoleculeModifierAndScorer)f;
		this.efunc = mof.getEfunc(d);
		this.dof = (FreeDihedral)mof.getDOFs().get(d);
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// pretend like each function call is a kernel launch
		// so they can only communicate via shared memory
		
		// init memory
		xdmin = mins.get(d);
		xdmax = maxs.get(d);
		xd = x.get(d);
		
		Residue res = dof.getResidue();
		coords = res.coords;
		dihedralAtomIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
		List<Integer> rotatedAtomIndicesList = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
		rotatedAtomIndices = new int[rotatedAtomIndicesList.size()];
		for (int i=0; i<rotatedAtomIndicesList.size(); i++) {
			rotatedAtomIndices[i] = rotatedAtomIndicesList.get(i);
		}
		
		// init first energy calc
		doEnergy = true;
		energyxd = xd;
		
		// INPUTS: xdmin, xdmax, xd, coords, dihedralAtomIndices, rotatedAtomIndices, doEnergy, energyxd
		
		pose();
		energy();
		search1();
		pose();
		energy();
		search2();
		pose();
		energy();
		search3();
		pose();
		energy();
		search4();
		pose();
		energy();
		search5();
		pose();
		energy();
		search6();
		for (int i=0; i<10; i++) {
			searchSurf();
			pose();
			energy();
		}
		search7();
		pose();
		energy();
		search8();
		pose();
		energy();
		search9();
		
		// OUTPUTS: step, xdstar
		
		// update the step
		lastStep = step;
		if (iteration == 0) {
			firstStep = lastStep;
		}

		// update the dofs and the conf
		x.set(d, xdstar);
		f.setDOF(d, xdstar);
		
		iteration++;
	}
	
	private void pose() {
		
		if (!doEnergy) {
			return;
		}
		
		// too easy... gotta make it harder =P
		//dof.apply(energyxd);
		
		// NOTE: simulate the vector types on the gpu to make porting the logic easier
		
		// get the old dihedral
		
		// get the four atom positions: a, b, c, d
		double[] a = make();
		double[] b = make();
		double[] c = make();
		double[] d = make();
		
		copy(coords, dihedralAtomIndices[0]*3, a);
		copy(coords, dihedralAtomIndices[1]*3, b);
		copy(coords, dihedralAtomIndices[2]*3, c);
		copy(coords, dihedralAtomIndices[3]*3, d);
		
		// translate so everything is centered on b
		subtractInPlace(a, b);
		subtractInPlace(c, b);
		subtractInPlace(d, b);
		
		// build a right orthnormal matrix [rx,ry,rz] where z is bc and ba points along x
		double[] rz = make(c);
		normalizeInPlace(rz);
		double[] rx = subtract(a, scale(c, dot(a, c)/dot(c, c)));
		normalizeInPlace(rx);
		double[] ry = cross(rz, rx);
		normalizeInPlace(ry);
		
		// use r^{-1} to rotate d into our axis-aligned space
		rotateInverse(d, rx, ry, rz);
		
		// look at the x,y coords of d to get the dihedral angle
		d[2] = 0;
		normalizeInPlace(d);
		double currentSin = d[1];
		double currentCos = d[0];
		
		// get the delta dihedral
		double newDihedral = Math.toRadians(energyxd);
		double newSin = Math.sin(newDihedral);
		double newCos = Math.cos(newDihedral);
		double deltaSin = newSin*currentCos - newCos*currentSin;
		double deltaCos = newCos*currentCos + newSin*currentSin;
		
		// build the delta rotation matrix about z by angle deltaDihedral
		double[] dx = { deltaCos, deltaSin, 0 };
		double[] dy = { -deltaSin, deltaCos, 0 };
		double[] dz = UnitZ;
		
		// apply transformation to every downstream atom
		double[] p = make();
		for (int index : rotatedAtomIndices) {
			copy(coords, index*3, p);
			subtractInPlace(p, b);
			rotateInverse(p, rx, ry, rz);
			rotate(p, dx, dy, dz);
			rotate(p, rx, ry, rz);
			addInPlace(p, b);
			copy(p, coords, index*3);
		}
	}
	
	private void energy() {
		
		if (!doEnergy) {
			return;
		}
		
		energyfxd = efunc.getEnergy();
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
