package edu.duke.cs.osprey.energy.forcefield;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.NEIGHBORTYPE;

public class BigForcefieldEnergy implements EnergyFunction {
	
	private static final long serialVersionUID = 5242606861996508290L;

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
	
	// physical constants and constant params
	private static final double coulombConstant = 332.0;
	private static final double solvCutoff = 9.0;
	
	private ForcefieldParams params;
	private Groups groups;
	
	// NOTE: use buffers here instead of arrays to make syncing with GPU easier
	
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
	
	// pre-computed vdW, electrostatics, and solvation params
	// layout per atom pair: Aij, Bij, charge, lambda1, radius1, alpha1, lambda2, radius2, alpha2
	private DoubleBuffer precomputed;
	
	private double internalSolvEnergy;
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions) {
		this(params, interactions, false);
	}
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions, boolean useDirectBuffers) {
		
		// TODO: implement solvation toggle
		// TODO: implement dynamic vs static atom groups
		
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
		if (useDirectBuffers) {
			coords = ByteBuffer.allocateDirect(numAtoms*3*Double.BYTES).order(ByteOrder.nativeOrder()).asDoubleBuffer();
		} else {
			coords = DoubleBuffer.allocate(numAtoms*3);
		}
		updateCoords();
		
		// do one pass over the group pairs to count the number of atom pairs
		int numAtomPairs = 0;
		for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
			
			AtomGroup[] groupPair = interactions.get(groupPairIndex);
			AtomGroup group1 = groupPair[0];
			AtomGroup group2 = groupPair[1];
			
			for (NEIGHBORTYPE type : Arrays.asList(NEIGHBORTYPE.BONDED14, NEIGHBORTYPE.NONBONDED)) {
				numAtomPairs += AtomNeighbors.getPairIndicesByType(
					group1.getAtoms(),
					group2.getAtoms(),
					group1 == group2,
					type
				).size();
			}
		}
		
		// NOTE: DAGK has 32,006,379 atom pairs
		// the precomputed buffer would take 2,304,459,288 bytes (2.15 GiB)
		// sadly, this actually overflows an int for the buffer size
		// overflows aren't well detected by java, so let's explicitly check here
		long bufferSize = (long)numAtomPairs * 9*Double.BYTES;
		if (bufferSize > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Too many atom pairs, can't allocate large enough buffer");
		}
		
		// allocate our buffers
		if (useDirectBuffers) {
			atomFlags = ByteBuffer.allocateDirect(numAtomPairs*2*Integer.BYTES).order(ByteOrder.nativeOrder()).asIntBuffer();
			precomputed = ByteBuffer.allocateDirect(numAtomPairs*9*Double.BYTES).order(ByteOrder.nativeOrder()).asDoubleBuffer();
		} else {
			atomFlags = IntBuffer.allocate(numAtomPairs*2);
			precomputed = DoubleBuffer.allocate(numAtomPairs*9);
		}
		num14Pairs = 0;
		numNbPairs = 0;
		
		NBParams nbparams1 = new NBParams();
		NBParams nbparams2 = new NBParams();
		VdwParams vdwparams = new VdwParams();
		SolvParams solvparams1 = new SolvParams();
		SolvParams solvparams2 = new SolvParams();
		
		// pre-pre-compute some vdW constants
		double Bmult = params.vdwMultiplier*params.vdwMultiplier;
		Bmult = Bmult*Bmult*Bmult;
		double Amult = Bmult*Bmult;
		
		// pre-pre-compute some solvation constants
		double solvCoeff = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI));
		solvCoeff *= params.solvScale;
		
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
					
					precomputed.put(vdwparams.Aij);
					precomputed.put(vdwparams.Bij);
					precomputed.put(atom1.charge*atom2.charge);
	
					// is this a heavy pair?
					if (!atom1.isHydrogen() && !atom2.isHydrogen()) {
						
						// save the solvation params
						getSolvParams(atom1, solvparams1);
						getSolvParams(atom2, solvparams2);
						
						double alpha1 = solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
						double alpha2 = solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
				
						precomputed.put(solvparams1.lambda);
						precomputed.put(solvparams1.radius);
						precomputed.put(alpha1);
						precomputed.put(solvparams2.lambda);
						precomputed.put(solvparams2.radius);
						precomputed.put(alpha2);
						
					} else {
						
						// use bogus info for the precalculated parameters
						// yeah, it takes up extra space, but space is cheap
						// it's more important that we make the location of these
						// parameters predictable so we can use parallelism to compute energies
						// also, syncing the solvation params with the vdw/es params helps cache performance
						for (int j=0; j<6; j++) {
							precomputed.put(0);
						}
					}
				}
			}
		}
		
		// flip our buffers
		atomFlags.flip();
		precomputed.flip();
		
		// compute internal solvation energy
		// ie, add up all the dGref terms for all atoms in internal groups
		internalSolvEnergy = 0;
		for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
			
			AtomGroup[] groupPair = interactions.get(groupPairIndex);
			AtomGroup group1 = groupPair[0];
			AtomGroup group2 = groupPair[1];
			
			if (group1 == group2) {
				for (Atom atom : group1.getAtoms()) {
					if (!atom.isHydrogen()) {
						getSolvParams(atom, solvparams1);
						internalSolvEnergy += solvparams1.dGref;
					}
				}
			}
		}
		internalSolvEnergy *= params.solvScale;
	}
	
	public DoubleBuffer getCoords() {
		return coords;
	}
	
	public IntBuffer getAtomFlags() {
		return atomFlags;
	}
	
	public DoubleBuffer getPrecomputed() {
		return precomputed;
	}
	
	public int getNumAtomPairs() {
		return num14Pairs + numNbPairs;
	}
	
	public int getNum14AtomPairs() {
		return num14Pairs;
	}
	
	public double getCoulombFactor() {
		return coulombConstant/params.dielectric;
	}
	
	public double getScaledCoulombFactor() {
		return getCoulombFactor()*params.forcefld.getCoulombScaling();
	}
	
	public double getSolvationCutoff2() {
		return solvCutoff*solvCutoff;
	}
	
	public boolean useDistDependentDielectric() {
		return params.distDepDielect;
	}
	
	public boolean useHElectrostatics() {
		return params.hElect;
	}
	
	public boolean useHVdw() {
		return params.hVDW;
	}
	
	public double getInternalSolvationEnergy() {
		return internalSolvEnergy;
	}
	
	public void updateCoords() {
		coords.rewind();
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			coords.put(group.getCoords());
		}
	}
	
	@Override
	public double getEnergy() {
		
		// OPTIMIZATION: this function gets hit a lot! so even pedantic optimizations can make a difference
		// I've also tweaked the code with fancy scoping to try to reduce register pressure
		// the idea is to limit the scope of temporary variables as much as possible
		// so the compiler/jvm has the most flexibilty to use registers
		// this actually has a measurable impact on performance
		
		// copy some things to the local stack
		int num14Pairs = this.num14Pairs;
		int numAtomPairs = getNumAtomPairs();
		boolean distDepDielect = useDistDependentDielectric();
		boolean useHEs = useHElectrostatics();
		boolean useHVdw = useHVdw();
		IntBuffer atomFlags = this.atomFlags;
		DoubleBuffer precomputed = this.precomputed;
		double coulombFactor = getCoulombFactor();
		double scaledCoulombFactor = getScaledCoulombFactor();
		double solvCutoff2 = getSolvationCutoff2();
		
		// the benchmarks seem to run faster when the loop-scoped variables are declared here
		boolean bothHeavy, inRangeForSolv;
		double r2;
		double r = 0;
		int i9;
		
		// compute all the energies
		double esEnergy = 0;
		double vdwEnergy = 0;
		double solvEnergy = internalSolvEnergy;
		atomFlags.rewind();
		precomputed.rewind();
		for (int i=0; i<numAtomPairs; i++) {
			i9 = i*9;
			
			// read flags
			{
				int atom1Index, atom2Index;
				{
					int i2 = i*2;
					int atom1Flags = atomFlags.get(i2);
					int atom2Flags = atomFlags.get(i2 + 1);
					atom1Index = getAtomIndex(atom1Flags);
					atom2Index = getAtomIndex(atom2Flags);
					bothHeavy = !isHydrogen(atom1Flags) && !isHydrogen(atom2Flags);
				}
				
				// get the squared radius
				{
					int atom1Index3 = atom1Index*3;
					int atom2Index3 = atom2Index*3;
					double rx = coords.get(atom1Index3) - coords.get(atom2Index3);
					double ry = coords.get(atom1Index3 + 1) - coords.get(atom2Index3 + 1);
					double rz = coords.get(atom1Index3 + 2) - coords.get(atom2Index3 + 2);
					r2 = rx*rx + ry*ry + rz*rz;
				}
			}
			
			// do we need the sqrt?
			// they're expensive to compute, so let's only do it once per atom pair
			inRangeForSolv = r2 < solvCutoff2;
			if (!distDepDielect || (bothHeavy && inRangeForSolv)) {
				r = Math.sqrt(r2);
			}
			
			if (bothHeavy || useHEs) {
				
				// read precomputed params
				double charge = precomputed.get(i9 + 2);
				
				// compute electrostatics
				boolean is14Pair = i < num14Pairs;
				esEnergy += (is14Pair ? scaledCoulombFactor : coulombFactor)
					/ (distDepDielect ? r2 : r)
					* charge;
			}

			if (bothHeavy || useHVdw) {
				
				// read precomputed params
				double Aij = precomputed.get(i9);
				double Bij = precomputed.get(i9 + 1);
				
				// compute vdw
				double r6 = r2*r2*r2;
				double r12 = r6*r6;
				vdwEnergy += Aij/r12 - Bij/r6;
			}
			
			if (bothHeavy && inRangeForSolv) {
					
				// read precomputed params
				double lambda1 = precomputed.get(i9 + 3);
				double radius1 = precomputed.get(i9 + 4);
				double alpha1 = precomputed.get(i9 + 5);
				double lambda2 = precomputed.get(i9 + 6);
				double radius2 = precomputed.get(i9 + 7);
				double alpha2 = precomputed.get(i9 + 8);
				
				// compute solvation energy
				double Xij = (r - radius1)/lambda1;
				double Xji = (r - radius2)/lambda2;
				solvEnergy -= (alpha1*Math.exp(-Xij*Xij) + alpha2*Math.exp(-Xji*Xji))/r2;
			}
		}
		
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
}
