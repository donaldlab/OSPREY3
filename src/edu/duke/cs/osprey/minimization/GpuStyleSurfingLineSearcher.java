package edu.duke.cs.osprey.minimization;

import static edu.duke.cs.osprey.tools.VectorAlgebra.*;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.List;

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.VectorAlgebra;

public class GpuStyleSurfingLineSearcher implements LineSearcher {
	
	// NOTE: this class isn't useful by itself
	// but it was a prototype for how the GPU-based version works
	// and it's probably easier to understand the logic by reading this one =)
	
	// lots of the methods are public so their results can be individually compared to the gpu implementation
	
	private static final double Tolerance = 1e-6;
	private static final double InitialStepSize = 0.25; // for dihedral dofs
	
	private ObjectiveFunction.OneDof f;
	private EnergyFunction efunc;
	
	private double firstStep;
	private double lastStep;
	private int iteration;
	
	private double xdmin;
	private double xdmax;
	private double xd;
	private double fxd;
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
	private boolean surfHitEdge;
	
	private double[] coords;
	private int[] dihedralAtomIndices;
	private int[] rotatedAtomIndices;
	
	private boolean doEnergy;
	private double energyxd;
	private double energyfxd;
	
	@Override
	public void init(ObjectiveFunction.OneDof f) {
		
		this.f = f;
		
		firstStep = 1;
		lastStep = 1;
		iteration = 0;
		
		// HACKHACK: break the interfaces to get the pieces we need
		// the runtime will complain if we can't do these casts, so we're still safe, but the compiler has no idea what's going on
		MoleculeModifierAndScorer mof = (MoleculeModifierAndScorer)f.getParent();
		this.efunc = mof.getEfunc(f.getDimension());
		FreeDihedral dof = (FreeDihedral)mof.getDOFs().get(f.getDimension());
		
		// get coords and atom indices
		Residue res = dof.getResidue();
		coords = res.coords;
		dihedralAtomIndices = res.template.getDihedralDefiningAtoms(dof.getDihedralNumber());
		List<Integer> rotatedAtomIndicesList = res.template.getDihedralRotatedAtoms(dof.getDihedralNumber());
		rotatedAtomIndices = new int[rotatedAtomIndicesList.size()];
		for (int i=0; i<rotatedAtomIndicesList.size(); i++) {
			rotatedAtomIndices[i] = rotatedAtomIndicesList.get(i);
		}
	}
	
	@Override
	public double search(double xd) {
		
		preSearch(xd);
		
		// pretend like each function call is a kernel launch
		// so they can only communicate via shared memory
		
		for (int i=0; i<=3; i++) {
			pose();
			energy();
			search(i);
		}
		for (int i=0; i<100; i++) {
			surf();
			if (i % 2 == 0 && !doEnergy) {
				break;
			}
			pose();
			energy();
		}
		search4();
		for (int i=5; i<=6; i++) {
			pose();
			energy();
			search(i);
		}
		pose();
		
		postSearch();
		
		return xdstar;
	}
	
	public void preSearch(double xd) {
		
		// zero out memory
		surfStep = 0;
		
		xdmin = 0;
		xdmax = 0;
		this.xd = 0;
		fxd = 0;
		step = 0;
		xdm = 0;
		xdp = 0;
		fxdm = 0;
		fxdp = 0;
		xdstar = 0;
		fxdstar = 0;
		stepScale = 0;
		xdsurfHere = 0;
		fxdsurfHere = 0;
		energyxd = 0;
		energyfxd = 0;
		
		isSurfing = false;
		surfHitEdge = false;
		
		// get the step size, try to make it adaptive (based on historical steps if possible; else on step #)
		if (Math.abs(lastStep) > Tolerance && Math.abs(firstStep) > Tolerance) {
			step = InitialStepSize*Math.abs(lastStep / firstStep);
		} else {
			step = InitialStepSize/Math.pow(iteration + 1, 3);
		}
		
		// init args
		this.xd = xd;
		xdmin = f.getXMin();
		xdmax = f.getXMax();
		doEnergy = true;
		energyxd = xd;
	}
	
	public void pose() {
		
		if (!doEnergy) {
			return;
		}
		
		// can't call this on the GPU, gotta do it ourselves
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
	
	public void energy() {
		
		if (!doEnergy) {
			return;
		}
		
		energyfxd = efunc.getEnergy();
	}
	
	public void search(int i) {
		
		// this is dumb, but whatever...
		// no function pointers ftw!
		switch (i) {
			case 0: search0(); break;
			case 1: search1(); break;
			case 2: search2(); break;
			case 3: search3(); break;
			case 4: search4(); break;
			case 5: search5(); break;
			case 6: search6(); break;
		}
	}
	
	public void search0() {
		
		// read energy result
		fxd = energyfxd;
		
		// get the minus neighbor
		xdm = xd - step;
		
		// init next energy calc
		doEnergy = xdm >= xdmin;
		energyxd = xdm;
	}
	
	public void search1() {
		
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
	
	public void search2() {
		
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
	
	public void search3() {
		
		// read energy result
		fxdstar = energyfxd;
		
		// if we went downhill, use growing steps
		if (fxdstar < fxd + getTolerance(fxd)) {
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
	
	public void surf() {
		
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
			
			// if we hit the edge, stop anyway
			if (surfHitEdge) {
				isSurfing = false;
			}
		}
			
		if (isSurfing) {
		
			// surf a step
			double xdsurfNext = xd + stepScale*(xdsurfHere - xd);
			
			// clamp the x val
			if (xdsurfNext < xdmin) {
				surfHitEdge = true;
				xdsurfNext = xdmin;
			} else if (xdsurfNext > xdmax) {
				surfHitEdge = true;
				xdsurfNext = xdmax;
			}
			
			// init energy calc for new surf spot
			surfStep++;
			doEnergy = true;
			energyxd = xdsurfNext;
		}
	}
	
	public void search4() {
		
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
	
	public void search5() {
		
		double xdstarOld = xdstar;
		
		// read energy result
		if (doEnergy) {
			if (energyfxd < fxdstar - getTolerance(fxdstar)) {
				xdstar = energyxd;
				fxdstar = energyfxd;
			}
		}
		
		// init next energy calc
		double xdp = xdstarOld + 1;
		doEnergy = xdp <= xdmax;
		energyxd = xdp;
	}
	
	public void search6() {
		
		// read energy result
		if (doEnergy) {
			if (energyfxd < fxdstar - getTolerance(fxdstar)) {
				xdstar = energyxd;
				fxdstar = energyfxd;
			}
		}
		
		// update the step to the one we actually took
		step = xdstar - xd;
		
		// prep one final pose to leave the protein in the final pose
		doEnergy = true;
		energyxd = xdstar;
	}
	
	public void postSearch() {
		
		// update the step
		lastStep = step;
		if (iteration == 0) {
			firstStep = lastStep;
		}
		
		iteration++;
	}
	
	private double getTolerance(double f) {
		
		// use full tolerance, unless f is very small
		// then scale by the magnitude of f
		return Tolerance * Math.max(1, Math.abs(f));
	}
	
	public void check(String name, DoubleBuffer coords, int coordsOffset, ByteBuffer args, double internalSolvEnergy, double firstStep, double lastStep, int iteration) {
		try {
			
			// check coords
			double[] cpuCoord = make();
			double[] gpuCoord = make();
			for (int index : rotatedAtomIndices) {
				copy(this.coords, index*3, cpuCoord);
				copy(coords, (index + coordsOffset)*3, gpuCoord);
				check(cpuCoord, gpuCoord);
			}
			
			// check fields
			check(this.firstStep, firstStep);
			check(this.lastStep, lastStep);
			check(this.iteration, iteration);
			
			// check args
			check(rotatedAtomIndices.length, args.getInt(0));
			check(surfStep, args.getInt(4));
			check(xdmin, args.getDouble(8));
			check(xdmax, args.getDouble(16));
			check(xd, args.getDouble(24));
			check(fxd, args.getDouble(32));
			check(step, args.getDouble(40));
			check(xdm, args.getDouble(48));
			check(xdp, args.getDouble(56));
			check(fxdm, args.getDouble(64));
			check(fxdp, args.getDouble(72));
			check(xdstar, args.getDouble(80));
			check(fxdstar, args.getDouble(88));
			check(stepScale, args.getDouble(96));
			check(xdsurfHere, args.getDouble(104));
			check(fxdsurfHere, args.getDouble(112));
			check(energyxd, args.getDouble(120));
			check(internalSolvEnergy, args.getDouble(128));
			check(isSurfing, args.get(136) == 1);
			check(surfHitEdge, args.get(137) == 1);
			
		} catch (Throwable t) {
			throw new Error("check failure after kernel: " + name, t);
		}
	}
	
	public void check(double exp, double obs) {
		check(exp, obs, 1e-4);
	}
	
	public void check(double exp, double obs, double epsilon) {
		double absErr = Math.abs(exp - obs);
		if (absErr > epsilon) {
			throw new Error(String.format("expected: %12.6f, observed: %12.6f, absErr: %.12f", exp, obs, absErr));
		}
	}
	
	public void check(double[] exp, double[] obs) {
		check(exp, obs, 1e-4);
	}
	
	public void check(double[] exp, double[] obs, double epsilon) {
		double dist = VectorAlgebra.distance(exp, obs);
		if (dist > epsilon) {
			throw new Error(String.format("\nexpected: (%12.6f,%12.6f,%12.6f)\nobserved: (%12.6f,%12.6f,%12.6f)\ndist:     %.12f",
				exp[0], exp[1], exp[2], obs[0], obs[1], obs[2], dist
			));
		}
	}
	
	public void check(int exp, int obs) {
		if (exp != obs) {
			throw new Error(String.format("expected: %d, observed: %d", exp, obs));
		}
	}
	
	public void check(boolean exp, boolean obs) {
		if (exp != obs) {
			throw new Error(String.format("expected: %b, observed: %b", exp, obs));
		}
	}
	
	public void compare(double exp, double obs) {
		System.out.println(String.format("%12.6f, %12.6f, diff: %.16f, equal: %b", exp, obs, exp - obs, exp == obs));
	}
}
