package edu.duke.cs.osprey.energy.forcefield;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;

public class BigForcefieldEnergy {
	
	private static class Groups {
		
		private List<AtomGroup> groups;
		private Map<Integer,Integer> groupIndicesById;
		private int[] groupIndicesByPair;
		
		public Groups(int numPairs) {
			groups = new ArrayList<>();
			groupIndicesById = new HashMap<>();
			groupIndicesByPair = new int[numPairs*2];
		}
		
		public void addPair(int pairIndex, AtomGroup group1, AtomGroup group2) {
			groupIndicesByPair[pairIndex*2 + 0] = addOrGetGroupIndex(group1);
			groupIndicesByPair[pairIndex*2 + 1] = addOrGetGroupIndex(group2);
		}
		
		private int addOrGetGroupIndex(AtomGroup group) {
			
			// do we have an index already?
			Integer groupIndex = groupIndicesById.get(group.getId());
			if (groupIndex != null) {
				return groupIndex;
			}
			
			// nope, make one
			groupIndex = groups.size();
			groups.add(group);
			groupIndicesById.put(group.getId(), groupIndex);
			
			return groupIndex;
		}
		
		public int getNumGroups() {
			return groups.size();
		}
		
		public AtomGroup get(int index) {
			return groups.get(index);
		}
		
		public int getGroup1Index(int pairIndex) {
			return groupIndicesByPair[pairIndex*2 + 0];
		}
		
		public int getGroup2Index(int pairIndex) {
			return groupIndicesByPair[pairIndex*2 + 1];
		}
	}
	
	private static class VdwParams {
		public double Aij;
		public double Bij;
	}
	
	private ForcefieldParams params;
	private ForcefieldInteractions interactions;
	private Groups groups;
	
	// atom coordinates for all groups
	private int[] atomOffsets;
	// layout per atom: x, y, z
	private DoubleBuffer coords;
	
	// pre-computed vdW parameters
	// layout per atom pair:
	//    ints:    atom1 flags, atom2 flags
	//    doubles: Aij, Bij
	private IntBuffer flags14;
	private DoubleBuffer pre14;
	private IntBuffer flagsNb;
	private DoubleBuffer preNb;
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions) {
		
		// make sure the settings are supported
		if (!params.hElect) {
			throw new UnsupportedOperationException("toggling eletrostatics for hydrogens is not yet supported by this forcefield implementation");
		}
		if (!params.hVDW) {
			throw new UnsupportedOperationException("toggling vdW for hydrogens is not yet supported by this forcefield implementation");
		}
		// TODO: implement these if anyone cares?
		
		// TODO: implement solvation precomputation
		// TODO: implement energy calculation!
		
		this.params = params;
		this.interactions = interactions;
		
		// get one list of the unique atom groups in a stable order
		// (this is all the variable info, collecting it in one place will make uploading to the gpu faster)
		groups = new Groups(interactions.size());
		for (int i=0; i<interactions.size(); i++) {
			groups.addPair(i, interactions.get(i)[0], interactions.get(i)[1]);
		}
		
		// convert the group list into an atom list
		int numAtoms = 0;
		atomOffsets = new int[groups.getNumGroups()];
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			atomOffsets[i] = numAtoms;
			numAtoms += group.getAtoms().size();
		}
		coords = DoubleBuffer.allocate(numAtoms*3);
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			coords.put(group.getCoords());
		}
		coords.flip();
		
		// pre-pre-compute some vdW constants
		double Bmult = params.vdwMultiplier*params.vdwMultiplier;
		Bmult = Bmult*Bmult*Bmult;
		double Amult = Bmult*Bmult;
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		VdwParams vdwparams = new VdwParams();
		
		// precompute all the position-independent params
		
		// first, vdW params
		for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
			
			AtomGroup[] groupPair = interactions.get(groupPairIndex);
			AtomGroup group1 = groupPair[0];
			AtomGroup group2 = groupPair[1];
			int group1Index = groups.getGroup1Index(groupPairIndex);
			int group2Index = groups.getGroup2Index(groupPairIndex);
			
			// handle 1-4 interactions
			List<int[]> atomPairs14 = AtomNeighbors.getPairIndicesByType(
				group1.getAtoms(),
				group2.getAtoms(),
				group1 == group2,
				AtomNeighbors.NEIGHBORTYPE.BONDED14
			);
			
			flags14 = IntBuffer.allocate(atomPairs14.size()*2);
			pre14 = DoubleBuffer.allocate(atomPairs14.size()*2);
			
			for (int i=0; i<atomPairs14.size(); i++) {
				int[] atomIndices = atomPairs14.get(i);
				Atom atom1 = group1.getAtoms().get(atomIndices[0]);
				Atom atom2 = group2.getAtoms().get(atomIndices[1]);
				
				// save the atom flags
				flags14.put(makeFlags(group1Index, atomIndices[0], atom1));
				flags14.put(makeFlags(group2Index, atomIndices[1], atom2));
				
				// save the vdw params
				getNonBondedParams(atom1, nbparams1);
				getNonBondedParams(atom2, nbparams2);
				calcVdw(nbparams1, nbparams2, Amult, Bmult, vdwparams);
				
				// vdW scaling for 1-4 interactions
				vdwparams.Aij *= params.forcefld.getAij14Factor();
				vdwparams.Bij *= params.forcefld.getBij14Factor();
				
				pre14.put(vdwparams.Aij);
				pre14.put(vdwparams.Bij);
			}
			
			// handle non-bonded interactions
			List<int[]> atomPairsNb = AtomNeighbors.getPairIndicesByType(
				group1.getAtoms(),
				group2.getAtoms(),
				group1 == group2,
				AtomNeighbors.NEIGHBORTYPE.NONBONDED
			);
			
			flagsNb = IntBuffer.allocate(atomPairsNb.size()*2);
			preNb = DoubleBuffer.allocate(atomPairsNb.size()*2);
			
			for (int i=0; i<atomPairsNb.size(); i++) {
				int[] atomIndices = atomPairsNb.get(i);
				Atom atom1 = group1.getAtoms().get(atomIndices[0]);
				Atom atom2 = group2.getAtoms().get(atomIndices[1]);
				
				// save the atom flags
				flagsNb.put(makeFlags(group1Index, atomIndices[0], atom1));
				flagsNb.put(makeFlags(group2Index, atomIndices[1], atom2));
				
				// save the vdw params
				getNonBondedParams(atom1, nbparams1);
				getNonBondedParams(atom2, nbparams2);
				calcVdw(nbparams1, nbparams2, Amult, Bmult, vdwparams);
				preNb.put(vdwparams.Aij);
				preNb.put(vdwparams.Bij);
			}
		}
		
		// second, solvation params
		
		// pseudocode:
		// for each heavy atom pair 1-4 bonded or farther
		//    save pairs of solvation param lists
		// the lists don't actually interact,
		// but we can duplicate some values to increase cpu cache performance
	}

	private int makeFlags(int groupIndex, int atomIndexInGroup, Atom atom) {
		return makeFlags(
			atomOffsets[groupIndex] + atomIndexInGroup,
			atom.elementType.equalsIgnoreCase("H")
		);
	}

	private int makeFlags(int atomIndex, boolean isHydrogen) {
		
		// we could use fancy bit-wise encoding if we need it,
		// but since we only have one int and one boolean,
		// this is waaaaay easier
		
		// atomIndex can be zero, which doesn't have a sign, so bump it up one
		if (atomIndex == Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Really?? We have billions of atoms??");
		}
		atomIndex++;
		
		if (isHydrogen) {
			return atomIndex;
		} else {
			return -atomIndex;
		}
	}
	
	public int getAtomIndex(int flags) {
		// undo the bump we did in makeFlags14()
		return Math.abs(flags) - 1;
	}
	
	public boolean getIsHydrogen(int flags) {
		assert (flags != 0);
		return flags > 0;
	}
	
	private void calcVdw(NBParams nbparams1, NBParams nbparams2, double Amult, double Bmult, VdwParams vdwparams) {
		
		// Aij = (ri+rj)^12 * sqrt(ei*ej)
		// Bij = (ri+rj)^6 * sqrt(ei*ej)
		
		double epsilon = Math.sqrt(nbparams1.epsilon*nbparams2.epsilon);
		double radiusSum = nbparams1.r + nbparams2.r;
		vdwparams.Bij = radiusSum*radiusSum;
		vdwparams.Bij = vdwparams.Bij*vdwparams.Bij*vdwparams.Bij;
		vdwparams.Aij = vdwparams.Bij*vdwparams.Bij;
		vdwparams.Aij *= epsilon*Amult;
		vdwparams.Bij *= epsilon*Bmult;
	}
	
	private void getNonBondedParams(Atom atom, NBParams nbparams) {
		
		// HACKHACK: overrides for C atoms
		if (atom.elementType.equalsIgnoreCase("C") && params.forcefld.reduceCRadii()) {
			
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
}
