package edu.duke.cs.osprey.energy.forcefield;

import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.Map;

import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.VdwParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairs;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;

public class ResPairCache {
	
	public static class ResPair {
		
		public final int resIndex1;
		public final int resIndex2;
		public final Residue res1;
		public final Residue res2;
		public final AtomPairInfo info;
		public final double solvEnergy;
		public final double weight;
		public final double offset;
		
		public ResPair(int resIndex1, int resIndex2, Residue res1, Residue res2, double weight, double offset, AtomPairInfo info, SolvationForcefield.ResiduesInfo solvInfo, double solvScale) {
			
			this.resIndex1 = resIndex1;
			this.resIndex2 = resIndex2;
			this.res1 = res1;
			this.res2 = res2;
			this.info = info;

			// update the res pair offset with solvation energy if needed
			if (solvInfo != null) {
				solvEnergy = solvInfo.getResPairEnergy(res1, res2)*solvScale;
			} else {
				solvEnergy = 0.0;
			}
			
			this.weight = weight;
			this.offset = offset;
		}
	}
	
	public static class AtomPairInfo {
		
		public final int numAtomPairs;
		
		// layout per atom pair: flags (bit is14Bonded, bit isHeavyPair, 6 bits space, 3 bytes space, short atom1Offset, short atom2Offset)
		// yeah, there's some extra space here, but keeping 8-byte alignments is fast on 64 bit machines
		public final long[] flags;
		
		// layout per atom pair: charge, Aij, Bij, radius1, lambda1, alpha1, radius2, lambda2, alpha2
		public final double[] precomputed;
		public final int numPrecomputedPerAtomPair;
		
		public AtomPairInfo(Residue res1, Residue res2, ForcefieldParams ffparams, AtomPairs atomPairs, SolvationForcefield.ResiduesInfo solvInfo) {
			
			VdwParams vdwparams = new VdwParams();
			
			// how many precomputed values per atom pair?
			int numPrecomputedPerAtomPair = 1 + 2;
			if (solvInfo != null) {
				numPrecomputedPerAtomPair += solvInfo.getNumPrecomputedPerAtomPair();
			}
			this.numPrecomputedPerAtomPair = numPrecomputedPerAtomPair;
			
			// build the atom pairs
			
			// count the number of atom pairs and allocate space 
			this.numAtomPairs = atomPairs.getNumPairs(AtomNeighbors.Type.BONDED14) + atomPairs.getNumPairs(AtomNeighbors.Type.NONBONDED);
			this.flags = new long[numAtomPairs];
			this.precomputed = new double[numAtomPairs*numPrecomputedPerAtomPair];
			
			int flagsIndex = 0;
			int precomputedIndex = 0;
			
			for (AtomNeighbors.Type type : Arrays.asList(AtomNeighbors.Type.BONDED14, AtomNeighbors.Type.NONBONDED)) {
				for (int[] atomPair : atomPairs.getPairs(type)) {
			
					Atom atom1 = res1.atoms.get(atomPair[0]);
					Atom atom2 = res2.atoms.get(atomPair[1]);
					
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
					this.flags[flagsIndex++] = flags;
					
					// calc electrostatics params
					precomputed[precomputedIndex++] = atom1.charge*atom2.charge;
					
					// calc vdw params
					ffparams.getVdwParams(atom1, atom2, type, vdwparams);
					precomputed[precomputedIndex++] = vdwparams.Aij;
					precomputed[precomputedIndex++] = vdwparams.Bij;
					
					// calc solvation params if needed
					if (solvInfo != null) {
						solvInfo.putPrecomputed(precomputed, precomputedIndex, res1, atomPair[0], res2, atomPair[1], ffparams.solvScale);
						precomputedIndex += solvInfo.getNumPrecomputedPerAtomPair();
					}
				}
			}
		}
	}
	
	public final ForcefieldParams ffparams;
	public final AtomConnectivity connectivity;
	
	private Map<AtomPairs,AtomPairInfo> infos;
	
	public ResPairCache(ForcefieldParams ffparams, AtomConnectivity connectivity) {
		this.ffparams = ffparams;
		this.connectivity = connectivity;
		this.infos = new IdentityHashMap<>();
	}
	
	public ResPair get(Residues residues, ResidueInteractions.Pair pair, SolvationForcefield.ResiduesInfo solvInfo) {
		
		// lookup the residues
		int indxe1 = residues.findIndex(pair.resNum1);
		int index2 = residues.findIndex(pair.resNum2);
		Residue res1 = residues.get(indxe1);
		Residue res2 = residues.get(index2);
		
		AtomPairs atomPairs = connectivity.getAtomPairs(res1, res2);
		if (atomPairs == null) {
			throw new RuntimeException("Atom connectivity was not correctly calculated."
					+ " Can't find atom pairs for residues: " + res1.fullName + ", " + res2.fullName);
		}
		
		// look in the cache
		AtomPairInfo info = infos.get(atomPairs);
		if (info == null) {
			
			// cache miss!
			info = new AtomPairInfo(
				res1, res2,
				ffparams,
				atomPairs,
				solvInfo
			);
			
			infos.put(atomPairs, info);
		}
	
		return new ResPair(
			indxe1, index2,
			res1, res2,
			pair.weight, pair.offset,
			info,
			solvInfo,
			ffparams.solvScale
		);
	}
}
