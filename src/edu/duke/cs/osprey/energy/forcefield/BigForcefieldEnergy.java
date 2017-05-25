package edu.duke.cs.osprey.energy.forcefield;

import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.EEF1.SolvParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions.AtomGroup;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.NBParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.structure.Atom;
import edu.duke.cs.osprey.structure.AtomNeighbors;
import edu.duke.cs.osprey.structure.AtomNeighbors.Type;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class BigForcefieldEnergy implements EnergyFunction.DecomposableByDof, EnergyFunction.ExplicitChemicalChanges {
	
	private static final long serialVersionUID = 5242606861996508290L;
	
	public static class ParamInfo {
		
		// physical constants and constant params
		public static boolean printWarnings = true;
		private static final double coulombConstant = 332.0;
		private static final double solvCutoff = 9.0;
		
		public final ForcefieldParams params;
		public final double Bmult;
		public final double Amult;
		public final double solvCoeff;
		public final double coulombFactor;
		public final double scaledCoulombFactor;
		public final double solvationCutoff2;
		public final boolean useDistDependentDielectric;
		public final boolean useHElectrostatics;
		public final boolean useHVdw;
		public final boolean useEEF1;
		
		public ParamInfo(ForcefieldParams params) {
			
			this.params = params;
			
			// pre-pre-compute some vdW constants
			double vdw2 = params.vdwMultiplier*params.vdwMultiplier;
			this.Bmult = vdw2*vdw2*vdw2;
			this.Amult = Bmult*Bmult;
			
			// pre-pre-compute some solvation constants
			this.solvCoeff = 2.0/(4.0*Math.PI*Math.sqrt(Math.PI)) * params.solvScale;
			
			// pre-pre-compute some more stuff
			this.coulombFactor = coulombConstant/params.dielectric;
			this.scaledCoulombFactor = this.coulombFactor*params.forcefld.coulombScaling;
			this.solvationCutoff2 = solvCutoff*solvCutoff;
			
			// save some flags
			this.useDistDependentDielectric = params.distDepDielect;
			this.useHElectrostatics = params.hElect;
			this.useHVdw = params.hVDW;
			this.useEEF1 = params.solvationForcefield == SolvationForcefield.EEF1;
		}
	}
	
	private static class GroupPair {
		
		private static class Entry {
			
			private static final int NumFlagsPerPair = 2;
			private static final int NumPrecomputedPerPair = 9;
			
			public final List<int[]> atomPairs;
			public final boolean[] flags;
			public final double[] precomputed;
			
			public int atomPairOffset;
			
			public Entry(List<int[]> atomPairs) {
				this.atomPairs = atomPairs;
				this.flags = new boolean[atomPairs.size()*NumFlagsPerPair];
				this.precomputed = new double[atomPairs.size()*NumPrecomputedPerPair];
				this.atomPairOffset = -1;
			}
			
			public void addFlags(int pairIndex, boolean isHydrogen1, boolean isHydrogen2) {
				int i = pairIndex*NumFlagsPerPair;
				flags[i++] = isHydrogen1;
				flags[i++] = isHydrogen2;
			}
			
			public void addPrecomputed(int pairIndex, double Aij, double Bij, double charge, double lambda1, double radius1, double alpha1, double lambda2, double radius2, double alpha2) {
				int i = pairIndex*NumPrecomputedPerPair;
				precomputed[i++] = Aij;
				precomputed[i++] = Bij;
				precomputed[i++] = charge;
				precomputed[i++] = lambda1;
				precomputed[i++] = radius1;
				precomputed[i++] = alpha1;
				precomputed[i++] = lambda2;
				precomputed[i++] = radius2;
				precomputed[i++] = alpha2;
			}
			
			public boolean isHydrogen1(int pairIndex) {
				return flags[pairIndex*2];
			}
			
			public boolean isHydrogen2(int pairIndex) {
				return flags[pairIndex*2 + 1];
			}
		}
		
		private AtomGroup group1;
		private AtomGroup group2;
		private int sequenceNumber1;
		private int sequenceNumber2;
		private Map<Type,Entry> entries;
		private double internalSolvEnergy;
		
		public GroupPair(ParamInfo pinfo, AtomGroup group1, AtomGroup group2) {
			this.group1 = group1;
			this.group2 = group2;
			this.entries = new EnumMap<>(Type.class);
			build(pinfo);
		}
		
		public Entry get(Type type) {
			return entries.get(type);
		}
		
		public int getNumAtomPairs() {
			int count = 0;
			for (Entry entry : entries.values()) {
				count += entry.atomPairs.size();
			}
			return count;
		}
		
		public double getInternalSolvationEnergy() {
			return internalSolvEnergy;
		}
		
		public boolean handleChemicalChanges(ParamInfo pinfo) {
			
			// are there any chemical changes?
			if (group1.getSequenceNumber() == sequenceNumber1 && group2.getSequenceNumber() == sequenceNumber2) {
				return false;
			}
			
			build(pinfo);
			return true;
		}
		
		public void build(ParamInfo pinfo) {
			
			// build the bond type entries
			entries.clear();
			for (Type type : Arrays.asList(Type.BONDED14, Type.NONBONDED)) {
				entries.put(type, build(pinfo, type));
			}
			
			internalSolvEnergy = 0;
			if (pinfo.useEEF1) {
				
				// calc the internal solv energy
				// ie, add up all the dGref terms for all atoms in internal groups
				SolvParams solvparams = new SolvParams();
					
				if (group1 == group2) {
					for (Atom atom : group1.getAtoms()) {
						if (!atom.isHydrogen()) {
							pinfo.params.eef1parms.getSolvationParametersOrDefaults(atom, solvparams);
							internalSolvEnergy += solvparams.dGref;
						}
					}
				}
				internalSolvEnergy *= pinfo.params.solvScale;
			}
			
			sequenceNumber1 = group1.getSequenceNumber();
			sequenceNumber2 = group2.getSequenceNumber();
		}
		
		private Entry build(ParamInfo pinfo, Type type) {
			
			List<int[]> atomPairs = AtomNeighbors.getPairIndicesByType(
				group1.getAtoms(),
				group2.getAtoms(),
				group1 == group2,
				type
			);
			
			NBParams nbparams1 = new NBParams();
			NBParams nbparams2 = new NBParams();
			VdwParams vdwparams = new VdwParams();
			SolvParams solvparams1 = new SolvParams();
			SolvParams solvparams2 = new SolvParams();
			
			Entry entry = new Entry(atomPairs);
			
			// precompute static vars for each atom pair
			for (int i=0; i<atomPairs.size(); i++) {
				int[] atomIndices = atomPairs.get(i);
				Atom atom1 = group1.getAtoms().get(atomIndices[0]);
				Atom atom2 = group2.getAtoms().get(atomIndices[1]);
				
				entry.addFlags(i, atom1.isHydrogen(), atom2.isHydrogen());
				
				// save the vdw params
				pinfo.params.getNonBondedParametersOrThrow(atom1, type, nbparams1);
				pinfo.params.getNonBondedParametersOrThrow(atom2, type, nbparams2);
				calcVdw(nbparams1, nbparams2, pinfo.Amult, pinfo.Bmult, vdwparams);
				
				// vdW scaling for 1-4 interactions
				if (type == Type.BONDED14) {
					vdwparams.Aij *= pinfo.params.forcefld.Aij14Factor;
					vdwparams.Bij *= pinfo.params.forcefld.Bij14Factor;
				} else if (type == Type.NONBONDED) {
					vdwparams.Bij *= 2;
				}
				
				// is this a heavy pair?
				if (!atom1.isHydrogen() && !atom2.isHydrogen()) {
					
					// save the solvation params
					pinfo.params.eef1parms.getSolvationParametersOrDefaults(atom1, solvparams1);
					pinfo.params.eef1parms.getSolvationParametersOrDefaults(atom2, solvparams2);
					
					double alpha1 = pinfo.solvCoeff*solvparams1.dGfree*solvparams2.volume/solvparams1.lambda;
					double alpha2 = pinfo.solvCoeff*solvparams2.dGfree*solvparams1.volume/solvparams2.lambda;
			
					entry.addPrecomputed(i,
						vdwparams.Aij, vdwparams.Bij,
						atom1.charge*atom2.charge,
						solvparams1.lambda, solvparams1.radius, alpha1,
						solvparams2.lambda, solvparams2.radius, alpha2
					);
					
				} else {
					
					// don't need to compute solvation for this pair, just pad with zeros
					entry.addPrecomputed(i,
						vdwparams.Aij, vdwparams.Bij,
						atom1.charge*atom2.charge,
						0, 0, 0,
						0, 0, 0
					);
				}
			}
			
			return entry;
		}
	}

	private class Groups {
		
		private List<AtomGroup> groups;
		private int sequenceNumber;
		private Map<Integer,Integer> groupIndicesById;
		private int[] groupIndicesByPair;
		private GroupPair[] pairs;
		private Map<Integer,Integer> pairIndicesByIds;
		
		public Groups(int numPairs) {
			groups = new ArrayList<>();
			sequenceNumber = 0;
			groupIndicesById = new HashMap<>();
			groupIndicesByPair = new int[numPairs*2];
			pairs = new GroupPair[numPairs];
			pairIndicesByIds = new HashMap<>();
		}
		
		public void addPair(int pairIndex, GroupPair pair) {
			groupIndicesByPair[pairIndex*2 + 0] = addOrGetGroupIndex(pair.group1);
			groupIndicesByPair[pairIndex*2 + 1] = addOrGetGroupIndex(pair.group2);
			pairs[pairIndex] = pair;
			boolean wasAdded = pairIndicesByIds.put(makePairKey(pair.group1, pair.group2), pairIndex) == null;
			if (!wasAdded) {
				throw new IllegalArgumentException("group pair was already added, can't add again");
			}
		}
		
		public int getNumGroupPairs() {
			return pairs.length;
		}
		
		public GroupPair getPair(int pairIndex) {
			return pairs[pairIndex];
		}
		
		public GroupPair getPair(AtomGroup[] groupPair) {
			int key = makePairKey(groupPair[0], groupPair[1]);
			Integer pairIndex = pairIndicesByIds.get(key);
			if (pairIndex == null) {
				throw new IllegalArgumentException("group pair not found in this forcefield");
			}
			return getPair(pairIndex);
		}
		
		private int makePairKey(AtomGroup group1, AtomGroup group2) {
			checkId(group1.getId());
			checkId(group2.getId());
			return (group1.getId() << 16) | group2.getId();
		}
		
		private void checkId(int id) {
			if (id < 0 | id > 32767) {
				throw new IllegalArgumentException("group id " + id + " out of range [0,32767]");
			}
		}
		
		public Integer getGroupIndex(AtomGroup group) {
			return groupIndicesById.get(group.getId());
		}
		
		private int addOrGetGroupIndex(AtomGroup group) {
			
			// do we have an index already?
			Integer groupIndex = getGroupIndex(group);
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
		
		public int handleChemicalChanges() {
			
			boolean hasChange = false;
			for (GroupPair pair : pairs) {
				hasChange |= pair.handleChemicalChanges(pinfo);
			}
			
			if (hasChange) {
				sequenceNumber++;
			}
			
			return sequenceNumber;
		}
		
		public int getGroup1Index(int pairIndex) {
			return groupIndicesByPair[pairIndex*2 + 0];
		}
		
		public int getGroup2Index(int pairIndex) {
			return groupIndicesByPair[pairIndex*2 + 1];
		}
		
		public List<GroupPair> match(ForcefieldInteractions interactions) {
			List<GroupPair> groupPairs = new ArrayList<>();
			for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
				groupPairs.add(getPair(interactions.get(groupPairIndex)));
			}
			return groupPairs;
		}
	}
	
	private static class VdwParams {
		public double Aij;
		public double Bij;
	}
	
	private ParamInfo pinfo;
	private ForcefieldInteractions interactions;
	private BufferTools.Type bufferType;
	private Groups groups;
	private int groupsSequenceNumber;
	
	// NOTE: use buffers here instead of arrays to make syncing with GPU easier
	
	// atom coordinates for all groups
	private int[] atomOffsets;
	// layout per atom: x, y, z
	private DoubleBuffer coords;
	
	// atom pair into
	// layout per atom pair: atom1 flags, atom2 flags
	private IntBuffer atomFlags;
	
	// pre-computed vdW, electrostatics, and solvation params
	// layout per atom pair: Aij, Bij, charge, lambda1, radius1, alpha1, lambda2, radius2, alpha2
	private DoubleBuffer precomputed;
	
	private Subset fullSubset;
	private Map<Residue,Subset> subsetCache;
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions) {
		this(params, interactions, BufferTools.Type.Normal);
	}
	
	public BigForcefieldEnergy(ForcefieldParams params, ForcefieldInteractions interactions, BufferTools.Type bufferType) {
		
		// TODO: implement solvation toggle
		
		this.pinfo = new ParamInfo(params);
		this.interactions = interactions;
		this.bufferType = bufferType;
		
		// compute all the info for each group pair
		groups = new Groups(interactions.size());
		for (int i=0; i<interactions.size(); i++) {
			groups.addPair(i, new GroupPair(pinfo, interactions.get(i)[0], interactions.get(i)[1]));
		}
		
		build();
		
		// make the full subset
		fullSubset = new Subset(interactions, false);
		
		subsetCache = new HashMap<>();
	}
	
	private void build() {
		
		// update the sequence number
		groupsSequenceNumber = groups.handleChemicalChanges();
		
		// convert the group list into an atom list
		int numAtoms = 0;
		atomOffsets = new int[groups.getNumGroups()];
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			atomOffsets[i] = numAtoms;
			numAtoms += group.getAtoms().size();
		}
		coords = makeOrResizeBuffer(coords, numAtoms*3);
		
		// do one pass over the group pairs to count the number of atom pairs
		int numAtomPairs = 0;
		for (int groupPairIndex=0; groupPairIndex<groups.getNumGroupPairs(); groupPairIndex++) {
			numAtomPairs += groups.getPair(groupPairIndex).getNumAtomPairs();
		}
		
		// NOTE: DAGK has 32,006,379 atom pairs
		// the precomputed buffer would take 2,304,459,288 bytes (2.15 GiB)
		// sadly, this actually overflows an int for the buffer size
		// overflows aren't well detected by java, so let's explicitly check here
		long bufferSize = (long)numAtomPairs*GroupPair.Entry.NumPrecomputedPerPair;
		if (bufferSize > Integer.MAX_VALUE) {
			throw new IllegalArgumentException("Too many atom pairs, can't allocate large enough buffer");
		}
		
		// allocate our buffers
		atomFlags = makeOrResizeBuffer(atomFlags, numAtomPairs*GroupPair.Entry.NumFlagsPerPair);
		precomputed = makeOrResizeBuffer(precomputed, numAtomPairs*GroupPair.Entry.NumPrecomputedPerPair);
		
		int atomPairOffset = 0;
		for (Type type : Arrays.asList(Type.BONDED14, Type.NONBONDED)) {
			
			// concatenate all the atom flags and precomputed values
			for (int groupPairIndex=0; groupPairIndex<interactions.size(); groupPairIndex++) {
				GroupPair groupPair = groups.getPair(groupPairIndex);
				GroupPair.Entry entry = groupPair.get(type);
				int group1Index = groups.getGroup1Index(groupPairIndex);
				int group2Index = groups.getGroup2Index(groupPairIndex);
			
				// update atom pair offset
				entry.atomPairOffset = atomPairOffset;
				atomPairOffset += entry.atomPairs.size();
				
				// append the flags
				for (int atomPairIndex=0; atomPairIndex<entry.atomPairs.size(); atomPairIndex++) {
					int[] atomIndices = entry.atomPairs.get(atomPairIndex);
					int atom1Index = getGlobalAtomIndex(group1Index, atomIndices[0]);
					int atom2Index = getGlobalAtomIndex(group2Index, atomIndices[1]);
					atomFlags.put(packFlags(atom1Index, entry.isHydrogen1(atomPairIndex)));
					atomFlags.put(packFlags(atom2Index, entry.isHydrogen2(atomPairIndex)));
				}
				
				// append the precomputed
				precomputed.put(entry.precomputed);
			}
		}
		
		assert (atomFlags.position() == numAtomPairs*GroupPair.Entry.NumFlagsPerPair);
		assert (precomputed.position() == numAtomPairs*GroupPair.Entry.NumPrecomputedPerPair);
		
		// flip our buffers
		atomFlags.flip();
		precomputed.flip();
	}
	
	private DoubleBuffer makeOrResizeBuffer(DoubleBuffer buf, int size) {
		if (buf == null || buf.capacity() < size) {
			buf = bufferType.makeDouble(size);
		} else {
			buf.clear();
		}
		assert (buf.capacity() >= size);
		return buf;
	}
	
	private IntBuffer makeOrResizeBuffer(IntBuffer buf, int size) {
		if (buf == null || buf.capacity() < size) {
			buf = bufferType.makeInt(size);
		} else {
			buf.clear();
		}
		assert (buf.capacity() >= size);
		return buf;
	}

	public ParamInfo getParams() {
		return pinfo;
	}
	
	public ForcefieldInteractions getInteractions() {
		return interactions;
	}
	
	public DoubleBuffer getCoords() {
		return coords;
	}
	
	public int getAtomOffset(Residue res) {
		return getAtomOffset(interactions.getResidueAtomGroup(res));
	}
	
	public int getAtomOffset(AtomGroup group) {
		Integer groupIndex = groups.getGroupIndex(group);
		if (groupIndex == null) {
			throw new IllegalArgumentException("group not found");
		}
		return getGlobalAtomIndex(groupIndex, 0);
	}
	
	public IntBuffer getAtomFlags() {
		return atomFlags;
	}
	
	public DoubleBuffer getPrecomputed() {
		return precomputed;
	}
	
	public Subset getFullSubset() {
		return fullSubset;
	}
	
	@Override
	public int handleChemicalChanges() {
		return fullSubset.handleChemicalChanges();
	}
	
	private int handleChemicalChangesInternal() {
		
		// if there was a chemical change, rebuild the forcefield
		int seq = groups.handleChemicalChanges();
		if (seq != groupsSequenceNumber) {
			groupsSequenceNumber = seq;
			build();
		}
		
		return seq;
	}
	
	public void updateCoords() {
		
		// rebuild if we have chemical changes
		handleChemicalChanges();
		
		// copy atom coords into the buffer
		coords.rewind();
		for (int i=0; i<groups.getNumGroups(); i++) {
			AtomGroup group = groups.get(i);
			coords.put(group.getCoords());
		}
	}
	
	@Override
	public double getEnergy() {
		return fullSubset.getEnergy();
	}
		
	@Override
	public List<EnergyFunction> decomposeByDof(Molecule m, List<DegreeOfFreedom> dofs) {
		
		List<EnergyFunction> efuncs = new ArrayList<>();
		
		for (DegreeOfFreedom dof : dofs) {

			Residue res = dof.getResidue();
			if (res == null) {
				
				// when there's no residue at the dof, then use the whole efunc
				efuncs.add(this);
				
			} else {
				
				// otherwise, make an efunc for only that residue
				// but share efuncs between dofs in the same residue
				Subset efunc = subsetCache.get(res);
				if (efunc == null) {
					efunc = new Subset(interactions.makeSubsetByResidue(res));
					subsetCache.put(res, efunc);
				}
				efuncs.add(efunc);
			}
		}
		
		return efuncs;
	}
	
	private static void calcVdw(NBParams nbparams1, NBParams nbparams2, double Amult, double Bmult, VdwParams vdwparams) {
	
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
	
	private int getGlobalAtomIndex(int groupIndex, int atomIndexInGroup) {
		return atomOffsets[groupIndex] + atomIndexInGroup;
	}

	private int packFlags(int atomIndex, boolean isHydrogen) {
		
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
	
	private int unpackAtomIndex(int flags) {
		assert (flags != 0);
		
		// undo the bump we did in makeFlags()
		return Math.abs(flags) - 1;
	}
	
	private boolean unpackIsHydrogen(int flags) {
		assert (flags != 0);
		return flags > 0;
	}
	
	public class Subset implements EnergyFunction.ExplicitChemicalChanges {
		
		private static final long serialVersionUID = -2038116438543332018L;
		
		private boolean makeTable;
		private List<GroupPair> groupPairs;
		private int sequenceNumber;
		private int numPairs;
		private int num14Pairs;
		private double internalSolvEnergy;
		private IntBuffer subsetTable;
		
		public Subset(ForcefieldInteractions interactions) {
			this(interactions, true);
		}
		
		public Subset(ForcefieldInteractions interactions, boolean makeTable) {
			
			this.makeTable = makeTable;
			
			// map interactions to group pairs
			groupPairs = groups.match(interactions);
			
			build();
		}
		
		private void build() {
			
			// count the atom pairs and the internal solvation energy
			numPairs = 0;
			num14Pairs = 0;
			internalSolvEnergy = 0;
			for (GroupPair pair : groupPairs) {
				numPairs += pair.getNumAtomPairs();
				num14Pairs += pair.get(Type.BONDED14).atomPairs.size();
				internalSolvEnergy += pair.getInternalSolvationEnergy();
			}
			
			if (makeTable) {
				
				// make the subset table
				subsetTable = makeOrResizeBuffer(subsetTable, numPairs);
				for (Type type : Arrays.asList(Type.BONDED14, Type.NONBONDED)) {
					for (GroupPair pair : groupPairs) {
						GroupPair.Entry entry = pair.get(type);
						for (int i=0; i<entry.atomPairs.size(); i++) {
							subsetTable.put(entry.atomPairOffset + i);
						}
					}
				}
				subsetTable.flip();
			}
		}
		
		public IntBuffer getSubsetTable() {
			return subsetTable;
		}
		
		public int getNumAtomPairs() {
			return numPairs;
		}
		
		public int getNum14AtomPairs() {
			return num14Pairs;
		}
		
		public double getInternalSolvationEnergy() {
			return internalSolvEnergy;
		}
		
		public boolean isBroken() {
			for (GroupPair pair : groupPairs) {
				if (pair.group1.isBroken() || pair.group2.isBroken()) {
					return true;
				}
			}
			return false;
		}
		
		@Override
		public int handleChemicalChanges() {
			int seq = handleChemicalChangesInternal();
			if (sequenceNumber != seq) {
				sequenceNumber = seq;
				build();
			}
			return seq;
		}
		
		@Override
		public double getEnergy() {
			
			// OPTIMIZATION: this function gets hit a lot! so even pedantic optimizations can make a difference
			// I've also tweaked the code with fancy scoping to try to reduce register pressure
			// the idea is to limit the scope of temporary variables as much as possible
			// so the compiler/jvm has the most flexibilty to use registers
			// this actually has a measurable impact on performance
			
			handleChemicalChanges();
			updateCoords();
			
			if (isBroken()) {
				return Double.POSITIVE_INFINITY;
			}
			
			// copy some things to the local stack
			IntBuffer subsetTable = this.subsetTable;
			int num14Pairs = this.num14Pairs;
			int numAtomPairs = this.numPairs;
			boolean distDepDielect = pinfo.useDistDependentDielectric;
			boolean useHEs = pinfo.useHElectrostatics;
			boolean useHVdw = pinfo.useHVdw;
			IntBuffer atomFlags = BigForcefieldEnergy.this.atomFlags;
			DoubleBuffer precomputed = BigForcefieldEnergy.this.precomputed;
			double coulombFactor = pinfo.coulombFactor;
			double scaledCoulombFactor = pinfo.scaledCoulombFactor;
			double solvCutoff2 = pinfo.solvationCutoff2;
			
			// the benchmarks seem to run faster when the loop-scoped variables are declared here
			boolean bothHeavy, inRangeForSolv;
			double r2;
			double r = 0;
			int i;
			int i2;
			int i9;
			
			// compute all the energies
			double esEnergy = 0;
			double vdwEnergy = 0;
			double solvEnergy = internalSolvEnergy;
			for (int j=0; j<numAtomPairs; j++) {
				
				if (subsetTable != null) {
					i = subsetTable.get(j);
				} else {
					i = j;
				}
				
				i2 = i*2;
				i9 = i*9;
				
				// read flags
				{
					int atom1Index, atom2Index;
					{
						int atom1Flags = atomFlags.get(i2);
						int atom2Flags = atomFlags.get(i2 + 1);
						atom1Index = unpackAtomIndex(atom1Flags);
						atom2Index = unpackAtomIndex(atom2Flags);
						bothHeavy = !unpackIsHydrogen(atom1Flags) && !unpackIsHydrogen(atom2Flags);
					}
					
					// get the squared radius
					{
						int atom1Index3 = atom1Index*3;
						int atom2Index3 = atom2Index*3;
						double d = coords.get(atom1Index3) - coords.get(atom2Index3);
						r2 = d*d;
						d = coords.get(atom1Index3 + 1) - coords.get(atom2Index3 + 1);
						r2 += d*d;
						d = coords.get(atom1Index3 + 2) - coords.get(atom2Index3 + 2);
						r2 += d*d;
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
				
				if (pinfo.useEEF1 && bothHeavy && inRangeForSolv) {
						
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
	}
}
