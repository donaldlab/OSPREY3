package edu.duke.cs.osprey.energy.forcefield;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.Type;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class ResidueForcefieldEnergy implements EnergyFunction.DecomposableByDof {
	
	private static final long serialVersionUID = -4768384219061898745L;
	
	private static final double coulombConstant = 332.0;
	private static final double solvCutoff = 9.0;
	private static final double solvCutoff2 = solvCutoff*solvCutoff;
	
	private static class VdwParams {
		public double Aij;
		public double Bij;
	}
	
	private ForcefieldParams params;
	private ResidueInteractions inters;
	private Molecule mol;
	
	private boolean isBroken;
	
	private double Bmult;
	private double Amult;
	private double solvCoeff;
	private double coulombFactor;
	private double scaledCoulombFactor;
	
	public ResidueForcefieldEnergy(ForcefieldParams params, ResidueInteractions inters, Molecule mol) {
		
		this.params = params;
		this.inters = inters;
		this.mol = mol;
		
		// is this a broken conformation?
		isBroken = false;
		for (ResidueInteractions.Pair pair : inters) {
			int numConfProblems = getRes(mol, pair.resNum1).confProblems.size()
				+ getRes(mol, pair.resNum2).confProblems.size();
			if (numConfProblems > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		
		// pre-pre-compute some constants
		double vdw2 = params.vdwMultiplier*params.vdwMultiplier;
		Bmult = vdw2*vdw2*vdw2;
		Amult = Bmult*Bmult;
		// TODO: make this a constant
		solvCoeff = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI)) * params.solvScale;
		coulombFactor = coulombConstant/params.dielectric;
		scaledCoulombFactor = coulombFactor*params.forcefld.coulombScaling;
	}
	
	private Residue getRes(Molecule mol, String resNum) {
		return mol.getResByPDBResNumber(resNum);
	}

	@Override
	public double getEnergy() {
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		VdwParams vdwparams = new VdwParams();
		SolvParams solvparams1 = new SolvParams();
		SolvParams solvparams2 = new SolvParams();
		
		double energy = 0;
		
		// for each residue pair...
		for (ResidueInteractions.Pair resPair : inters) {
			Residue res1 = getRes(mol, resPair.resNum1);
			Residue res2 = getRes(mol, resPair.resNum2);
			
			double pairEnergy = 0;
		
			for (Type atomPairType : Arrays.asList(AtomNeighbors.Type.BONDED14, AtomNeighbors.Type.NONBONDED)) {
				
				// for each atom pair...
				for (int[] atomIndices : AtomNeighbors.getPairIndicesByType(res1.atoms, res2.atoms, res1 == res2, atomPairType)) {
					Atom atom1 = res1.atoms.get(atomIndices[0]);
					Atom atom2 = res2.atoms.get(atomIndices[1]);
					
					boolean isHeavyPair = !atom1.isHydrogen() && !atom2.isHydrogen();
	
					// get the squared radius
					double r2 = 0;
					{
						int atom1Index3 = atomIndices[0]*3;
						int atom2Index3 = atomIndices[1]*3;
						double d = res1.coords[atom1Index3] - res2.coords[atom2Index3];
						r2 = d*d;
						d = res1.coords[atom1Index3 + 1] - res2.coords[atom2Index3 + 1];
						r2 += d*d;
						d = res1.coords[atom1Index3 + 2] - res2.coords[atom2Index3 + 2];
						r2 += d*d;
					}
					
					double r = Math.sqrt(r2);
					
					// electrostatics
					if (isHeavyPair || params.hElect) {
						
						double coulomb = atomPairType == Type.BONDED14 ? scaledCoulombFactor : coulombFactor;
						double effectiveR = params.distDepDielect ? r2 : r;
						
						pairEnergy += coulomb*atom1.charge*atom2.charge/effectiveR;
					}
					
					// van der Waals
					if (isHeavyPair || params.hVDW) {
						
						getNonBondedParams(atom1, nbparams1);
						getNonBondedParams(atom2, nbparams2);
						calcVdw(nbparams1, nbparams2, vdwparams, atomPairType);
						
						// compute vdw
						double r6 = r2*r2*r2;
						double r12 = r6*r6;
						pairEnergy += vdwparams.Aij/r12 - vdwparams.Bij/r6;
					}

					// solvation (EEF1)
					if (isHeavyPair && params.solvationForcefield == SolvationForcefield.EEF1 && r2 < solvCutoff2) {
						
						// save the solvation params
						getSolvParams(atom1, solvparams1);
						getSolvParams(atom2, solvparams2);
						
						double alpha1 = solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
						double alpha2 = solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
				
						// compute solvation energy
						double Xij = (r - solvparams1.radius)/solvparams1.lambda;
						double Xji = (r - solvparams2.radius)/solvparams2.lambda;
						pairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
					}
				}
			}
			
			// add the internal solvation energy if needed
			if (res1 == res2 && params.solvationForcefield == SolvationForcefield.EEF1) {

				double internalSolvEnergy = 0;
				
				// add up all the dGref terms for all the atoms
				for (Atom atom : res1.atoms) {
					if (!atom.isHydrogen()) {
						getSolvParams(atom, solvparams1);
						internalSolvEnergy += solvparams1.dGref;
					}
				}
				
				pairEnergy += internalSolvEnergy*params.solvScale;
			}
			
			// apply weights and offsets
			energy += pairEnergy*resPair.weight;
			energy += resPair.offset;
		}
		
		return energy;
	}
	
	private void getNonBondedParams(Atom atom, NBParams nbparams) {
		
		// HACKHACK: overrides for C atoms
		if (atom.isCarbon() && params.forcefld.reduceCRadii) {
			
			// Jeff: shouldn't these settings be in a config file somewhere?
			nbparams.epsilon = 0.1;
			nbparams.r = 1.9;
			
		} else {
			
			boolean success = params.getNonBondedParameters(atom.type, nbparams);
			if (!success) {
				// TODO: what's the right error-handling behavior here?
				// skip any atom pairs without params and keep computing?
				// use default values for nbparams?
				// or crash and tell the user to fix the problem?
				throw new Error("couldn't find non-bonded parameters for atom type: " + atom.forceFieldType);
			}
		}
	}
	
	private void calcVdw(NBParams nbparams1, NBParams nbparams2, VdwParams vdwparams, Type atomPairType) {
	
		// Aij = (ri+rj)^12 * sqrt(ei*ej)
		// Bij = (ri+rj)^6 * sqrt(ei*ej)
		
		double epsilon = Math.sqrt(nbparams1.epsilon*nbparams2.epsilon);
		double radiusSum = nbparams1.r + nbparams2.r;
		vdwparams.Bij = radiusSum*radiusSum;
		vdwparams.Bij = vdwparams.Bij*vdwparams.Bij*vdwparams.Bij;
		vdwparams.Aij = vdwparams.Bij*vdwparams.Bij;
		vdwparams.Aij *= epsilon*Amult;
		vdwparams.Bij *= epsilon*Bmult;
		
		// vdW scaling for 1-4 interactions
		if (atomPairType == Type.BONDED14) {
			vdwparams.Aij *= params.forcefld.Aij14Factor;
			vdwparams.Bij *= params.forcefld.Bij14Factor;
		} else if (atomPairType == Type.NONBONDED) {
			vdwparams.Bij *= 2;
		}
	}
	
	private static Set<String> warnedAtomTypes = new HashSet<>();
	
	private void getSolvParams(Atom atom, SolvParams solvparams) {
		boolean success = params.eef1parms.getSolvationParameters(atom, solvparams);
		if (!success) {
			
			// if there's no params, don't crash, use defaults instead
			if (warnedAtomTypes.add(atom.forceFieldType)) {
				System.err.println("WARNING: couldn't find solvation parameters for atom type: " + atom.forceFieldType + ", using default values");
			}
			
			solvparams.dGref = 0;
			solvparams.dGfree = 0;
			solvparams.volume = 0;
			solvparams.lambda = 1;
			solvparams.radius = 0;
		}
	}

	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		// TODO Auto-generated method stub
		return null;
	}
}
