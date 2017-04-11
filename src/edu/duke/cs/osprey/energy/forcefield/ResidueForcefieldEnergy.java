package edu.duke.cs.osprey.energy.forcefield;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairs;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class ResidueForcefieldEnergy implements EnergyFunction.DecomposableByDof {
	
	private static final long serialVersionUID = -4768384219061898745L;
	
	private static final double coulombConstant = 332.0;
	private static final double solvCutoff = 9.0;
	private static final double solvCutoff2 = solvCutoff*solvCutoff;
	private static final double solvTrig = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI));
	
	public static class ResPair {
		
		public final Residue res1;
		public final Residue res2;
		
		public double weight = 0;
		public double offset = 0;
		public int numAtomPairs = 0;
		
		// layout per atom pair: flags (bit is14Bonded, bit isHeavyPair, 6 bits space, 3 bytes space, short atom1Offset, short atom2Offset)
		// yeah, there's some extra space here, but keeping 8-byte alignments is fast on 64 bit machines
		public long[] flags = null;
		
		// layout per atom pair: charge, Aij, Bij, radius1, lambda1, alpha1, radius2, lambda2, alpha2
		public double[] precomputed = null;
		
		public ResPair(Molecule mol, ResidueInteractions.Pair pair) {
			
			this.res1 = mol.getResByPDBResNumber(pair.resNum1);
			this.res2 = mol.getResByPDBResNumber(pair.resNum2);
			
			this.weight = pair.weight;
			this.offset = pair.offset;
		}
	}
	
	public final ForcefieldParams params;
	public final ResidueInteractions inters;
	public final Molecule mol;
	public final AtomConnectivity connectivity;
	
	private Residue[] residues;
	private ResPair[] resPairs;
	
	private boolean isBroken;
	
	private double coulombFactor;
	private double scaledCoulombFactor;
	
	/* TODO: move to gpu area
	private ByteBuffer coords;
	
	 * buffer layout:
	 * NOTE: try to use 8-byte alignments to be happy on 64-bit machines
	 * 
	 * for each residue pair:
	 *    long numAtomPairs
	 *    
	 *    for each atom pair:
	 *       short atomIndex1
	 *       short atomIndex2
	 *       int flags  (isHeavyPair, is14Bonded)
	 *       double charge
	 *       double Aij
	 *       double Bij
	 *       double radius1
	 *       double lambda1
	 *       double alpha1
	 *       double radius2
	 *       double lambda2
	 *       double alpha2
	 *       
	 *    double weight
	 *    double offset
	private ByteBuffer precomputed;
	
	private static final int ResPairBytes = Long.BYTES + Double.BYTES*2;
	private static final int AtomPairBytes = Integer.BYTES + Short.BYTES*2 + Double.BYTES*9;
	*/
	
	public ResidueForcefieldEnergy(ForcefieldParams params, ResidueInteractions inters, Molecule mol, AtomConnectivity connectivity) {
		
		this.params = params;
		this.inters = inters;
		this.mol = mol;
		this.connectivity = connectivity;
		
		int index;
		
		// map the residues and pairs to the molecule
		residues = new Residue[inters.getResidueNumbers().size()];
		index = 0;
		for (String resNum : inters.getResidueNumbers()) {
			residues[index++] = mol.getResByPDBResNumber(resNum);
		}
		resPairs = new ResPair[inters.size()];
		index = 0;
		for (ResidueInteractions.Pair pair : inters) {
			resPairs[index++] = new ResPair(mol, pair);
		}
		
		// is this a broken conformation?
		isBroken = false;
		for (ResPair pair : resPairs) {
			if (pair.res1.confProblems.size() + pair.res2.confProblems.size() > 0) {
				isBroken = true;
				
				// we're done here, no need to analyze broken conformations
				return;
			}
		}
		
		// pre-compute some constants needed by getEnergy()
		coulombFactor = coulombConstant/params.dielectric;
		scaledCoulombFactor = coulombFactor*params.forcefld.coulombScaling;
		
		// TODO: can make lookup table by template for this
		// pre-compute internal solvation energies if needed
		double[] internalSolvEnergies;
		switch (params.solvationForcefield) {
			
			case EEF1:
				internalSolvEnergies = new double[mol.residues.size()];
				Arrays.fill(internalSolvEnergies, Double.NaN);
				SolvParams solvparams = new SolvParams();
				for (Residue res : residues) {
					
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
		
		/* TODO: move into GPU area
		// count atoms and offsets for each residue
		int[] atomOffsetsByResIndex = new int[mol.residues.size()];
		Arrays.fill(atomOffsetsByResIndex, -1);
		int atomOffset = 0;
		int numAtoms = 0;
		for (String resNum : inters.getResidueNumbers()) {
			Residue res = residuesByNum.get(resNum);
			atomOffsetsByResIndex[res.indexInMolecule] = atomOffset;
			atomOffset += 3*Double.BYTES*res.atoms.size();
			numAtoms += res.atoms.size();
		}
		
		// make the coords buffer
		coords = bufferType.make(numAtoms*3*Double.BYTES);
		*/
		
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		SolvParams solvparams1 = new SolvParams();
		SolvParams solvparams2 = new SolvParams();
		
		double vdw2 = params.vdwMultiplier*params.vdwMultiplier;
		double Bmult = vdw2*vdw2*vdw2;
		double Amult = Bmult*Bmult;
		double solvCoeff = solvTrig*params.solvScale;
		
		// build the atom pairs
		for (ResPair pair : resPairs) {
			
			AtomPairs atomPairs = connectivity.getAtomPairs(mol, pair.res1, pair.res2);
			
			// count the number of atom pairs and allocate space 
			pair.numAtomPairs = atomPairs.getNumPairs(AtomNeighbors.Type.BONDED14) + atomPairs.getNumPairs(AtomNeighbors.Type.NONBONDED);
			pair.flags = new long[pair.numAtomPairs];
			pair.precomputed = new double[pair.numAtomPairs*9];
			
			int flagsIndex = 0;
			int precomputedIndex = 0;
			
			for (AtomNeighbors.Type type : Arrays.asList(AtomNeighbors.Type.BONDED14, AtomNeighbors.Type.NONBONDED)) {
				for (int[] atomPair : atomPairs.getPairs(type)) {
			
					Atom atom1 = pair.res1.atoms.get(atomPair[0]);
					Atom atom2 = pair.res2.atoms.get(atomPair[1]);
					
					// pack the flags
					boolean isHeavyPair = !atom1.isHydrogen() && !atom2.isHydrogen();
					boolean is14Bonded = type == AtomNeighbors.Type.BONDED14;
					int atomOffset1 = atomPair[0]*3;
					int atomOffset2 = atomPair[1]*3;
					long flags = is14Bonded ? 1 : 0;
					flags <<= 1;
					flags |= isHeavyPair ? 1 : 0;
					flags <<= 46;
					flags |= (atomOffset1) & 0xffffL;
					flags <<= 16;
					flags |= (atomOffset2) & 0xffffL;
					pair.flags[flagsIndex++] = flags;
					
					// calc electrostatics params
					pair.precomputed[precomputedIndex++] = atom1.charge*atom2.charge;
					
					// calc vdW params
					// Aij = (ri+rj)^12 * sqrt(ei*ej)
					// Bij = (ri+rj)^6 * sqrt(ei*ej)
					
					getNonBondedParams(atom1, nbparams1);
					getNonBondedParams(atom2, nbparams2);
					double epsilon = Math.sqrt(nbparams1.epsilon*nbparams2.epsilon);
					double radiusSum = nbparams1.r + nbparams2.r;
					double Bij = radiusSum*radiusSum;
					Bij = Bij*Bij*Bij;
					double Aij = Bij*Bij;
					Aij *= epsilon*Amult;
					Bij *= epsilon*Bmult;
					
					// vdW scaling by connectivity
					if (is14Bonded) {
						Aij *= params.forcefld.Aij14Factor;
						Bij *= params.forcefld.Bij14Factor;
					} else {
						Bij *= 2;
					}
					
					pair.precomputed[precomputedIndex++] = Aij;
					pair.precomputed[precomputedIndex++] = Bij;
					
					// compute solvation params if needed
					if (isHeavyPair) {
						
						getSolvParams(atom1, solvparams1);
						getSolvParams(atom2, solvparams2);
						double alpha1 = solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
						double alpha2 = solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
						
						pair.precomputed[precomputedIndex++] = solvparams1.radius;
						pair.precomputed[precomputedIndex++] = solvparams1.lambda;
						pair.precomputed[precomputedIndex++] = alpha1;
						pair.precomputed[precomputedIndex++] = solvparams2.radius;
						pair.precomputed[precomputedIndex++] = solvparams2.lambda;
						pair.precomputed[precomputedIndex++] = alpha2;
						
					} else {
						
						precomputedIndex += 6;
					}
				}
			}
			
			// update the pair offset with internal solvation energy if needed
			if (pair.res1 == pair.res2) {
				switch (params.solvationForcefield) {
					
					case EEF1:
						pair.offset += internalSolvEnergies[pair.res1.indexInMolecule];
					break;
					
					default:
						// do nothing
				}
			}
		}
	}
	
	private void getNonBondedParams(Atom atom, NBParams nbparams) {
		
		// HACKHACK: overrides for C atoms
		if (atom.isCarbon() && params.forcefld.reduceCRadii) {
			
			// TODO: move this into Forcefield, and do the same for other forcefield energy classes too
			// Jeff: shouldn't these settings be in a config file somewhere?
			nbparams.epsilon = 0.1;
			nbparams.r = 1.9;
			
		} else {
			
			boolean success = params.getNonBondedParameters(atom.type, nbparams);
			if (!success) {
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
	
	/* TODO: move to GPU area
	public void copyCoords() {
		coords.clear();
		DoubleBuffer doubleCoords = coords.asDoubleBuffer();
		for (String resNum : inters.getResidueNumbers()) {
			doubleCoords.put(residuesByNum.get(resNum).coords);
		}
		coords.clear();
	}
	*/

	@Override
	public double getEnergy() {
		return getEnergy(resPairs);
	}
	
	private double getEnergy(ResPair[] resPairs) {
		
		// check broken-ness first. easy peasy
		if (isBroken) {
			return Double.POSITIVE_INFINITY;
		}
		
		// copy stuff to the stack/registers, to improve CPU cache performance
		boolean useHEs = params.hElect;
		boolean useHvdW = params.hVDW;
		double coulombFactor = this.coulombFactor;
		double scaledCoulombFactor = this.scaledCoulombFactor;
		boolean distDepDielect = params.distDepDielect;
		boolean useEEF1 = params.solvationForcefield == SolvationForcefield.EEF1;
		
		double energy = 0;
		
		for (int i=0; i<resPairs.length; i++) {
			ResPair pair = resPairs[i];
			
			// copy pair values/references to the stack/registers
			// so we don't have to touch the pair object while inside the atom pair loop
			double[] coords1 = pair.res1.coords;
			double[] coords2 = pair.res2.coords;
			int numAtomPairs = pair.numAtomPairs;
			long[] flags = pair.flags;
			double[] precomputed = pair.precomputed;
			
			double resPairEnergy = 0;
			
			// for each atom pair...
			int j9 = 0;
			for (int j=0; j<numAtomPairs; j++) {
				
				// read the flags
				// NOTE: this is efficient, but destructive to the val
				long atomPairFlags = flags[j];
				int atomOffset2 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 16;
				int atomOffset1 = (int)(atomPairFlags & 0xffff);
				atomPairFlags >>= 46;
				boolean isHeavyPair = (atomPairFlags & 0x1) == 0x1;
				atomPairFlags >>= 1;
				boolean is14Bonded = (atomPairFlags & 0x1) == 0x1;
				
				// get the radius
				double r2;
				double r;
				{
					// read atom coords
					double x1 = coords1[atomOffset1];
					double y1 = coords1[atomOffset1 + 1];
					double z1 = coords1[atomOffset1 + 2];
					double x2 = coords2[atomOffset2];
					double y2 = coords2[atomOffset2 + 1];
					double z2 = coords2[atomOffset2 + 2];
					
					double d;
					
					d = x1 - x2;
					r2 = d*d;
					d = y1 - y2;
					r2 += d*d;
					d = z1 - z2;
					r2 += d*d;
					r = Math.sqrt(r2);
				}
				
				// read the bit flags
				
				// electrostatics
				if (isHeavyPair || useHEs) {
					double charge = precomputed[j9++];
					if (is14Bonded) {
						if (distDepDielect) {
							resPairEnergy += scaledCoulombFactor*charge/r2;
						} else {
							resPairEnergy += scaledCoulombFactor*charge/r;
						}
					} else {
						if (distDepDielect) {
							resPairEnergy += coulombFactor*charge/r2;
						} else {
							resPairEnergy += coulombFactor*charge/r;
						}
					}
				} else {
					j9++;
				}
				
				// van der Waals
				if (isHeavyPair || useHvdW) {
					
					double Aij = precomputed[j9++];
					double Bij = precomputed[j9++];
					
					// compute vdw
					double r6 = r2*r2*r2;
					double r12 = r6*r6;
					resPairEnergy += Aij/r12 - Bij/r6;
					
				} else {
					j9 += 2;
				}
				
				// solvation
				if (useEEF1 && isHeavyPair && r2 < solvCutoff2) {
							
					double radius1 = precomputed[j9++];
					double lambda1 = precomputed[j9++];
					double alpha1 = precomputed[j9++];
					double radius2 = precomputed[j9++];
					double lambda2 = precomputed[j9++];
					double alpha2 = precomputed[j9++];
					
					// compute solvation energy
					double Xij = (r - radius1)/lambda1;
					double Xji = (r - radius2)/lambda2;
					resPairEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
					
				} else {
					j9 += 6;
				}
			}
			
			// apply weights and offsets
			energy += resPairEnergy*pair.weight;
			energy += pair.offset;
		}
		
		return energy;
	}
	
	private ResPair[] makeResPairsSubset(Residue res) {
	
		// pass 1: count
		int num = 0;
		for (ResPair resPair : resPairs) {
			if (resPair.res1 == res || resPair.res2 == res) {
				num++;
			}
		}
		
		// pass 2: collect
		ResPair[] pairs = new ResPair[num];
		num = 0;
		for (ResPair resPair : resPairs) {
			if (resPair.res1 == res || resPair.res2 == res) {
				pairs[num++] = resPair;
			}
		}
		return pairs;
	}
	
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule mol, List<DegreeOfFreedom> dofs) {
		
		class Subset implements EnergyFunction {
			
			private static final long serialVersionUID = 4664215035458391734L;
			
			private ResPair[] resPairs;
			
			@Override
			public double getEnergy() {
				return ResidueForcefieldEnergy.this.getEnergy(resPairs);
			}
		}
		
		Map<Residue,Subset> cache = new HashMap<>();
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		for (DegreeOfFreedom dof : dofs) {
			Residue res = dof.getResidue();
			
			if (res == null) {
				
				// no res, just use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// make a subset energy function
				Subset subset = cache.get(res);
				if (subset == null) {
					subset = new Subset();
					subset.resPairs = makeResPairsSubset(res);
					cache.put(res, subset);
				}
				efuncs.add(subset);
			}
		}
		
		return efuncs;
	}
}
