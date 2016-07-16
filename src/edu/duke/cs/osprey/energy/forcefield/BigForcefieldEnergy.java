package edu.duke.cs.osprey.energy.forcefield;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.NEIGHBORTYPE;

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
	private Groups groups;
	
	// atom coordinates for all groups
	private int[] atomOffsets;
	// layout per atom: x, y, z
	private DoubleBuffer coords;
	private int numAtoms;
	
	// atom pair into
	// layout per atom pair: atom1 flags, atom2 flags
	private IntBuffer atomFlags;
	private int num14Pairs;
	private int numNbPairs;
	
	// pre-computed vdW, electrostatics params
	// layout per atom pair: Aij, Bij
	private DoubleBuffer preVdwEs;
	
	// pre-computed solvation params
	// layout per atom: lambda, radius, alpha
	private DoubleBuffer preSolv;
	private double internalSolvEnergy;
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions) {
		
		// make sure the settings are supported
		if (!params.hElect) {
			throw new UnsupportedOperationException("toggling eletrostatics for hydrogens is not yet supported by this forcefield implementation");
		}
		if (!params.hVDW) {
			throw new UnsupportedOperationException("toggling vdW for hydrogens is not yet supported by this forcefield implementation");
		}
		// TODO: implement hydrogen toggles (if anyone cares?)
		// TODO: implement solvation toggle
		// TODO: implement dynamic vs static atom groups
		// TODO: internal solvation energy?
		
		this.params = params;
		
		// get one list of the unique atom groups in a stable order
		// (this is all the variable info, collecting it in one place will make uploading to the gpu faster)
		groups = new Groups(interactions.size());
		for (int i=0; i<interactions.size(); i++) {
			groups.addPair(i, interactions.get(i)[0], interactions.get(i)[1]);
		}
		
		// convert the group list into an atom list
		numAtoms = 0;
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
		
		// do one pass over the group pairs to count the number of atom pairs
		int numAtomPairs = 0;
		int numHeavyAtomPairs = 0;
		for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
			
			AtomGroup[] groupPair = interactions.get(groupPairIndex);
			AtomGroup group1 = groupPair[0];
			AtomGroup group2 = groupPair[1];
			
			for (NEIGHBORTYPE type : Arrays.asList(NEIGHBORTYPE.BONDED14, NEIGHBORTYPE.NONBONDED)) {
			
				List<int[]> atomPairs = AtomNeighbors.getPairIndicesByType(
					group1.getAtoms(),
					group2.getAtoms(),
					group1 == group2,
					type
				);
				
				for (int i=0; i<atomPairs.size(); i++) {
					int[] atomIndices = atomPairs.get(i);
					Atom atom1 = group1.getAtoms().get(atomIndices[0]);
					Atom atom2 = group2.getAtoms().get(atomIndices[1]);
					
					numAtomPairs++;
					
					if (!atom1.isHydrogen() && !atom2.isHydrogen()) {
						numHeavyAtomPairs++;
					}
				}
			}
		}
		
		// pre-pre-compute some vdW constants
		double Bmult = params.vdwMultiplier*params.vdwMultiplier;
		Bmult = Bmult*Bmult*Bmult;
		double Amult = Bmult*Bmult;
		
		// pre-pre-compute some solvation constants
		double solvCoeff = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI));
		
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		VdwParams vdwparams = new VdwParams();
		SolvParams solvparams1 = new SolvParams();
		SolvParams solvparams2 = new SolvParams();
		
		// allocate our buffers
		atomFlags = IntBuffer.allocate(numAtomPairs*2);
		preVdwEs = DoubleBuffer.allocate(numAtomPairs*3);
		num14Pairs = 0;
		numNbPairs = 0;
		preSolv = DoubleBuffer.allocate(numHeavyAtomPairs*6);
		
		for (NEIGHBORTYPE type : Arrays.asList(NEIGHBORTYPE.BONDED14, NEIGHBORTYPE.NONBONDED)) {
			
			// do another pass over the groups to precompute all the position-independent params
			for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
				
				AtomGroup[] groupPair = interactions.get(groupPairIndex);
				AtomGroup group1 = groupPair[0];
				AtomGroup group2 = groupPair[1];
				int group1Index = groups.getGroup1Index(groupPairIndex);
				int group2Index = groups.getGroup2Index(groupPairIndex);
			
				List<int[]> atomPairs = AtomNeighbors.getPairIndicesByType(
					group1.getAtoms(),
					group2.getAtoms(),
					group1 == group2,
					type
				);
			
				for (int i=0; i<atomPairs.size(); i++) {
					int[] atomIndices = atomPairs.get(i);
					Atom atom1 = group1.getAtoms().get(atomIndices[0]);
					Atom atom2 = group2.getAtoms().get(atomIndices[1]);
					
					// save atom flags
					atomFlags.put(makeFlags(group1Index, atomIndices[0], atom1));
					atomFlags.put(makeFlags(group2Index, atomIndices[1], atom2));
					
					if (type == NEIGHBORTYPE.BONDED14) {
						num14Pairs++;
					} else if (type == NEIGHBORTYPE.NONBONDED) {
						numNbPairs++;
					}
					
					// save the vdw params
					getNonBondedParams(atom1, nbparams1);
					getNonBondedParams(atom2, nbparams2);
					calcVdw(nbparams1, nbparams2, Amult, Bmult, vdwparams);
					
					// vdW scaling for 1-4 interactions
					if (type == NEIGHBORTYPE.BONDED14) {
						vdwparams.Aij *= params.forcefld.getAij14Factor();
						vdwparams.Bij *= params.forcefld.getBij14Factor();
					} else if (type == NEIGHBORTYPE.NONBONDED) {
						vdwparams.Bij *= 2;
					}
					
					preVdwEs.put(vdwparams.Aij);
					preVdwEs.put(vdwparams.Bij);
					preVdwEs.put(atom1.charge*atom2.charge);
	
					// is this a heavy pair?
					if (!atom1.isHydrogen() && !atom2.isHydrogen()) {
						
						// save the solvation params
						getSolvParams(atom1, solvparams1);
						getSolvParams(atom2, solvparams2);
						
						double alpha1 = solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
						double alpha2 = solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
				
						preSolv.put(solvparams1.lambda);
						preSolv.put(solvparams1.radius);
						preSolv.put(alpha1);
						preSolv.put(solvparams2.lambda);
						preSolv.put(solvparams2.radius);
						preSolv.put(alpha2);
					}
				}
			}
		}
		
		// flip our buffers
		atomFlags.flip();
		preVdwEs.flip();
		preSolv.flip();
		
		// compute internal solvation energy
		// ie, add up all the dGref terms for all atoms
		internalSolvEnergy = 0;
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			for (Atom atom : group.getAtoms()) {
				if (!atom.isHydrogen()) {
					getSolvParams(atom, solvparams1);
					internalSolvEnergy += solvparams1.dGref;
				}
			}
		}
	}
	
	public int getNumAtomPairs() {
		return num14Pairs + numNbPairs;
	}
	
	public double calculateTotalEnergy() {
		
		// physical constants and constant params
		final double coulombConstant = 332.0;
		final double solvCutoff = 9.0;
		
		// pre-compute some more constants
		double coulombFactor = coulombConstant/params.dielectric;
		double scaledCoulombFactor = coulombFactor*params.forcefld.getCoulombScaling();
		double solvCutoff2 = solvCutoff*solvCutoff;
		
		// compute all the energies
		double esEnergy = 0;
		double vdwEnergy = 0;
		double solvEnergy = internalSolvEnergy;
		atomFlags.rewind();
		preVdwEs.rewind();
		preSolv.rewind();
		int numVdwEsPairs = num14Pairs + numNbPairs;
		for (int i=0; i<numVdwEsPairs; i++) {
			boolean is14Pair = i < num14Pairs;
			
			// read flags
			int atom1Flags = atomFlags.get();
			int atom2Flags = atomFlags.get();
			int atom1Index = getAtomIndex(atom1Flags);
			int atom2Index = getAtomIndex(atom2Flags);
			boolean atom1isH = isHydrogen(atom1Flags);
			boolean atom2isH = isHydrogen(atom2Flags);
			
			// read precomputed vdw/es params
			double Aij = preVdwEs.get();
			double Bij = preVdwEs.get();
			double charge = preVdwEs.get();
			
			// get the squared radius
			double r2 = calcr2(atom1Index, atom2Index);
			
			// TODO: implement hydrogen flags
			
			// compute electrostatics
			double c;
			if (is14Pair) {
				c = scaledCoulombFactor;
			} else {
				c = coulombFactor;
			}
			if (params.distDepDielect) {
				c /= r2;
			} else {
				c /= Math.sqrt(r2);
			}
			c *= charge;
			esEnergy += c;

			// compute vdw
			double r6 = r2*r2*r2;
			double r12 = r6*r6;
			vdwEnergy += Aij/r12 - Bij/r6;
			
			if (!atom1isH && !atom2isH) {
				
				// read precomputed solvation params
				double lambda1 = preSolv.get();
				double radius1 = preSolv.get();
				double alpha1 = preSolv.get();
				double lambda2 = preSolv.get();
				double radius2 = preSolv.get();
				double alpha2 = preSolv.get();
				
				if (r2 < solvCutoff2) {
					
					// compute solvation energy
					double r = Math.sqrt(r2);
					double Xij = (r - radius1)/lambda1;
					double Xji = (r - radius2)/lambda2;
					solvEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
				}
			}
		}
		
		// just in case...
		assert (atomFlags.position() == atomFlags.limit());
		assert (preVdwEs.position() == preVdwEs.limit());
		assert (preSolv.position() == preSolv.limit());
		
		solvEnergy *= params.solvScale;
		
		return esEnergy + vdwEnergy + solvEnergy;
	}

	private int makeFlags(int groupIndex, int atomIndexInGroup, Atom atom) {
		return makeFlags(
			atomOffsets[groupIndex] + atomIndexInGroup,
			atom.isHydrogen()
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
		// undo the bump we did in makeFlags()
		return Math.abs(flags) - 1;
	}
	
	public boolean isHydrogen(int flags) {
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
		if (atom.isCarbon() && params.forcefld.reduceCRadii()) {
			
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
	
	private void getSolvParams(Atom atom, SolvParams solvparams) {
		boolean success = params.eef1parms.getSolvationParameters(atom, solvparams);
		if (!success) {
			throw new Error("couldn't find solvation parameters for atom type: " + atom.forceFieldType);
		}
	}
	
	private double calcr2(int atom1Index, int atom2Index) {
		int atom1Index3 = atom1Index*3;
		int atom2Index3 = atom2Index*3;
		double rx = coords.get(atom1Index3) - coords.get(atom2Index3);
		double ry = coords.get(atom1Index3 + 1) - coords.get(atom2Index3 + 1);
		double rz = coords.get(atom1Index3 + 2) - coords.get(atom2Index3 + 2);
		return rx*rx + ry*ry + rz*rz;
	}
}
