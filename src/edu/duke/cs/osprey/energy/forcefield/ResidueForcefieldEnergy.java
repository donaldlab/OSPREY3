package edu.duke.cs.osprey.energy.forcefield;

import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPair;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairList;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.Type;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class ResidueForcefieldEnergy implements EnergyFunction.DecomposableByDof {
	
	private static final long serialVersionUID = -4768384219061898745L;
	
	private static final double coulombConstant = 332.0;
	private static final double solvCutoff = 9.0;
	private static final double solvCutoff2 = solvCutoff*solvCutoff;
	private static final double solvTrig = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI));
	
	public final ForcefieldParams params;
	public final ResidueInteractions inters;
	public final Molecule mol;
	public final AtomConnectivity connectivity;
	
	private boolean isBroken;
	
	// TODO: do we need these here?
	private double coulombFactor;
	private double scaledCoulombFactor;
	
	private double[] internalSolvEnergies;
	
	private ByteBuffer coords;
	
	/* buffer layout:
	 * 
	 * for each residue pair:
	 *    long numAtomPairs
	 *    double weight
	 *    double offset
	 *    
	 *    for each atom pair:
	 *       bool isHeavyPair
	 *       bool is14Bonded
	 *       6 bytes space
	 *       int atomsIndex1
	 *       int atomsIndex2
	 *       double charge
	 *       double Aij
	 *       double Bij
	 *       double radius1
	 *       double lambda1
	 *       double alpha1
	 *       double radius2
	 *       double lambda2
	 *       double alpha2
	 */
	private ByteBuffer precomputed;
	
	private static final int ResPairBytes = Long.BYTES + Double.BYTES*2;
	private static final int AtomPairBytes = Long.BYTES + Integer.BYTES*2 + Double.BYTES*9;
	
	public ResidueForcefieldEnergy(ForcefieldParams params, ResidueInteractions inters, Molecule mol, AtomConnectivity connectivity) {
		this(params, inters, mol, connectivity, BufferTools.Type.Normal);
	}
	
	public ResidueForcefieldEnergy(ForcefieldParams params, ResidueInteractions inters, Molecule mol, AtomConnectivity connectivity, BufferTools.Type bufferType) {
		
		this.params = params;
		this.inters = inters;
		this.mol = mol;
		this.connectivity = connectivity;
		
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
		
		// pre-compute some constants needed by getEnergy()
		coulombFactor = coulombConstant/params.dielectric;
		scaledCoulombFactor = coulombFactor*params.forcefld.coulombScaling;
		
		// pre-compute internal solvation energies if needed
		switch (params.solvationForcefield) {
			
			case EEF1:
				internalSolvEnergies = new double[mol.residues.size()];
				Arrays.fill(internalSolvEnergies, Double.NaN);
				SolvParams solvparams = new SolvParams();
				for (String resNum : inters.getResidueNumbers()) {
					Residue res = mol.getResByPDBResNumber(resNum);
					
					// add up all the dGref terms for all the atoms
					double energy = 0;
					for (Atom atom : res.atoms) {
						if (!atom.isHydrogen()) {
							getSolvParams(atom, solvparams);
							energy += solvparams.dGref;
						}
					}
					energy *= params.solvScale;
					
					internalSolvEnergies[res.indexInMolecule] = energy;
				}
			break;
			
			default:
				internalSolvEnergies = null;
		}
		
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[mol.residues.size()];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (String resNum : inters.getResidueNumbers()) {
			Residue res = getRes(mol, resNum);
			atomOffsetsByResIndex[res.indexInMolecule] = atomOffset;
			atomOffset += 3*Double.BYTES*res.atoms.size();
			numAtoms += res.atoms.size();
		}
		
		// make the coords buffer
		coords = bufferType.make(numAtoms*3*Double.BYTES);
		
		// pass 1: count the atom pairs
		long[] numAtomPairs = new long[inters.size()];
		int index = 0;
		int bufferBytes = 0;
		for (ResidueInteractions.Pair resPair : inters) {
			Residue res1 = getRes(mol, resPair.resNum1);
			Residue res2 = getRes(mol, resPair.resNum2);
			
			bufferBytes += ResPairBytes;
			
			// count the 1-4 bonded and non-bonded pairs
			long count = 0;
			for (AtomPair pair : connectivity.getAtomPairs(mol, res1, res2).pairs) {
				if (pair.type == AtomNeighbors.Type.BONDED14 || pair.type == AtomNeighbors.Type.NONBONDED) {
					count++;
					bufferBytes += AtomPairBytes;
				}
			}
			
			numAtomPairs[index++] = count;
		}
		
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		SolvParams solvparams1 = new SolvParams();
		SolvParams solvparams2 = new SolvParams();
		
		double vdw2 = params.vdwMultiplier*params.vdwMultiplier;
		double Bmult = vdw2*vdw2*vdw2;
		double Amult = Bmult*Bmult;
		double solvCoeff = solvTrig*params.solvScale;
		
		// pass 2: precompute everything
		precomputed = bufferType.make(bufferBytes);
		index = 0;
		for (ResidueInteractions.Pair resPair : inters) {
			Residue res1 = getRes(mol, resPair.resNum1);
			Residue res2 = getRes(mol, resPair.resNum2);
			
			// compute the residue pair offset
			double offset = resPair.offset;
			if (res1 == res2) {
				switch (params.solvationForcefield) {
					
					case EEF1:
						offset += internalSolvEnergies[res1.indexInMolecule];
					break;
					
					default:
						// do nothing
				}
			}
			
			// just in case
			assert (Double.isFinite(resPair.weight));
			assert (Double.isFinite(offset));
			
			// collect the residue pair precomputed values
			precomputed.putLong(numAtomPairs[index++]);
			precomputed.putDouble(resPair.weight);
			precomputed.putDouble(offset);
			
			AtomPairList atomPairs = connectivity.getAtomPairs(mol, res1, res2);
			for (int i=0; i<atomPairs.size(); i++) {
		
				// we only care about 1-4 bonded and non-bonded pairs
				AtomNeighbors.Type atomPairType = atomPairs.getType(i);
				if (atomPairType != AtomNeighbors.Type.BONDED14 && atomPairType != AtomNeighbors.Type.NONBONDED) {
					continue;
				}
				
				int atomIndex1 = atomPairs.getIndex1(res1, res2, i);
				int atomIndex2 = atomPairs.getIndex2(res1, res2, i);
				Atom atom1 = res1.atoms.get(atomIndex1);
				Atom atom2 = res2.atoms.get(atomIndex2);
				
				// compute flags
				boolean isHeavyPair = !atom1.isHydrogen() && !atom2.isHydrogen();
				boolean is14Bonded = atomPairType == AtomNeighbors.Type.BONDED14;
				long flags = 0
					| (isHeavyPair ? 0x1 : 0)
					| (is14Bonded ? 0x2 : 0);
				
				precomputed.putLong(flags);
				
				// get atoms indices
				precomputed.putInt(atomOffsetsByResIndex[res1.indexInMolecule] + atomIndex1*3*Double.BYTES);
				precomputed.putInt(atomOffsetsByResIndex[res2.indexInMolecule] + atomIndex2*3*Double.BYTES);
				
				// compute physical values
				getNonBondedParams(atom1, nbparams1);
				getNonBondedParams(atom2, nbparams2);
				
				// calc vdW params
				// Aij = (ri+rj)^12 * sqrt(ei*ej)
				// Bij = (ri+rj)^6 * sqrt(ei*ej)
				
				double epsilon = Math.sqrt(nbparams1.epsilon*nbparams2.epsilon);
				double radiusSum = nbparams1.r + nbparams2.r;
				double Bij = radiusSum*radiusSum;
				Bij = Bij*Bij*Bij;
				double Aij = Bij*Bij;
				Aij *= epsilon*Amult;
				Bij *= epsilon*Bmult;
				
				// vdW scaling by connectivity
				if (atomPairType == Type.BONDED14) {
					Aij *= params.forcefld.Aij14Factor;
					Bij *= params.forcefld.Bij14Factor;
				} else if (atomPairType == Type.NONBONDED) {
					Bij *= 2;
				}
				
				precomputed.putDouble(atom1.charge*atom2.charge);
				precomputed.putDouble(Aij);
				precomputed.putDouble(Bij);
				
				// compute solvation if needed
				if (isHeavyPair) {
					
					getSolvParams(atom1, solvparams1);
					getSolvParams(atom2, solvparams2);
					double alpha1 = solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
					double alpha2 = solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
					
					precomputed.putDouble(solvparams1.radius);
					precomputed.putDouble(solvparams1.lambda);
					precomputed.putDouble(alpha1);
					precomputed.putDouble(solvparams2.radius);
					precomputed.putDouble(solvparams2.lambda);
					precomputed.putDouble(alpha2);
					
				} else {
					
					precomputed.putDouble(0);
					precomputed.putDouble(0);
					precomputed.putDouble(0);
					precomputed.putDouble(0);
					precomputed.putDouble(0);
					precomputed.putDouble(0);
				}
			}
		}
		precomputed.flip();
	}
	
	private Residue getRes(Molecule mol, String resNum) {
		return mol.getResByPDBResNumber(resNum);
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
	
	public void copyCoords() {
		coords.clear();
		DoubleBuffer doubleCoords = coords.asDoubleBuffer();
		for (String resNum : inters.getResidueNumbers()) {
			doubleCoords.put(getRes(mol, resNum).coords);
		}
		coords.clear();
	}

	@Override
	public double getEnergy() {
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		copyCoords();
		
		double energy = 0;
		
		precomputed.rewind();
		for (int i=0; i<inters.size(); i++) {
			
			// read the res pair info
			long numAtomPairs = precomputed.getLong();
			double weight = precomputed.getDouble();
			double offset = precomputed.getDouble();
			
			double resPairEnergy = 0;
			
			for (int j=0; j<numAtomPairs; j++) {
				
				// read the flags
				long flags = precomputed.getLong();
				boolean isHeavyPair = (flags & 0x1) == 0x1;
				boolean is14Bonded = (flags & 0x2) == 0x2;
				
				// read atoms offsets
				int atomsOffset1 = precomputed.getInt();
				int atomsOffset2 = precomputed.getInt();
				
				// read the physical values
				double charge = precomputed.getDouble();
				double Aij = precomputed.getDouble();
				double Bij = precomputed.getDouble();
				double radius1 = precomputed.getDouble();
				double lambda1 = precomputed.getDouble();
				double alpha1 = precomputed.getDouble();
				double radius2 = precomputed.getDouble();
				double lambda2 = precomputed.getDouble();
				double alpha2 = precomputed.getDouble();
			
				// get the squared radius
				double r2;
				{
					coords.position(atomsOffset1);
					double x1 = coords.getDouble();
					double y1 = coords.getDouble();
					double z1 = coords.getDouble();
					
					coords.position(atomsOffset2);
					double x2 = coords.getDouble();
					double y2 = coords.getDouble();
					double z2 = coords.getDouble();
					
					double d;
					
					d = x1 - x2;
					r2 = d*d;
					d = y1 - y2;
					r2 += d*d;
					d = z1 - z2;
					r2 += d*d;
				}
				
				double r = Math.sqrt(r2);
				
				// electrostatics
				if (isHeavyPair || params.hElect) {
					
					double coulomb = is14Bonded ? scaledCoulombFactor : coulombFactor;
					double effectiveR = params.distDepDielect ? r2 : r;
					
					resPairEnergy += coulomb*charge/effectiveR;
				}
				
				// van der Waals
				if (isHeavyPair || params.hVDW) {
					
					// compute vdw
					double r6 = r2*r2*r2;
					double r12 = r6*r6;
					resPairEnergy += Aij/r12 - Bij/r6;
				}

				// solvation
				if (isHeavyPair && r2 < solvCutoff2) {
					
					switch (params.solvationForcefield) {
						
						case EEF1:
							
							// compute solvation energy
							double Xij = (r - radius1)/lambda1;
							double Xji = (r - radius2)/lambda2;
							resPairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
							
						break;
						
						default:
							// do nothing
					}
				}
			}
			
			// apply weights and offsets
			energy += resPairEnergy*weight;
			energy += offset;
		}
		
		return energy;
	}
	
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		// TODO Auto-generated method stub
		return null;
	}
}
