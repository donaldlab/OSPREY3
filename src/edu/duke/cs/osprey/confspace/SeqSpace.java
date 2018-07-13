package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.tools.MathTools;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;


/**
 * analogous to the conf space, but just for sequences
 */
public class SeqSpace implements Serializable {

	private static final long serialVersionUID = -7309869148482636562L;

	/**
	 * Combine multiple sequence spaces into a single sequence space
	 *
	 * If any position appears in multiple spaces, it must have the same residue types in all spaces
	 */
	public static SeqSpace union(List<SeqSpace> seqSpaces) {
		return seqSpaces.stream()
			.reduce((a, b) -> union(a, b))
			.orElseThrow(() -> new IllegalArgumentException("conf spaces list is empty"));
	}

	/**
	 * Combine multiple sequence spaces into a single sequence space
	 *
	 * If any residue appears in multiple spaces, it must have the same wild type and residue types in all spaces
	 */
	public static SeqSpace union(SeqSpace a, SeqSpace b) {

		SeqSpace u = new SeqSpace();

		// add all positions from a
		for (Position apos : a.positions) {
			u.makePos(
				apos.resNum,
				apos.wildType != null ? apos.wildType.name : null,
				apos.resTypes.stream()
					.map(resType -> resType.name)
					.collect(Collectors.toList())
			);
		}

		// add unique positions from b
		for (Position bpos : b.positions) {

			// if a position already exists at this residue, make sure it's the same sequence space
			Position upos = u.getPosition(bpos.resNum);
			if (upos != null) {

				// check the wild type
				String uWildType = upos.wildType != null ? upos.wildType.name : null;
				String bWildType = bpos.wildType != null ? bpos.wildType.name : null;
				if (!Objects.equals(uWildType, bWildType)) {
					throw new IllegalArgumentException(String.format(
						"the two positions at residue %s have different wild types: %s != %s",
						upos.resNum, uWildType, bWildType
					));
				}

				// check the res types
				Set<String> utypes = bpos.resTypes.stream()
					.map(resType -> resType.name)
					.collect(Collectors.toSet());
				Set<String> btypes = bpos.resTypes.stream()
					.map(resType -> resType.name)
					.collect(Collectors.toSet());
				if (!utypes.equals(btypes)) {
					throw new IllegalArgumentException(String.format(
						"the two positions at residue %s have different residue types: %s != %s",
						upos.resNum, utypes, btypes
					));
				}

				// the positions match, no need to add a new pos to u
				continue;
			}

			// all is well, add the pos
			u.makePos(
				bpos.resNum,
				bpos.wildType.name,
				bpos.resTypes.stream()
					.map(resType -> resType.name)
					.collect(Collectors.toList())
			);
		}

		return u;
	}

	public class Position implements Comparable<Position>, Serializable {

		private static final long serialVersionUID = -8317027786742083135L;

		public final SeqSpace seqSpace = SeqSpace.this;

		public final int index;
		public final String resNum;
		public final ResType wildType;
		public final List<ResType> resTypes;
		public final List<ResType> mutations;

		private final Map<String,ResType> resTypesByName;

		private Position(int index, String resNum, String wildType, List<String> resTypes) {

			this.index = index;
			this.resNum = resNum;

			// make the res types
			this.resTypes = new ArrayList<>(resTypes.size());
			for (String resType : resTypes) {
				this.resTypes.add(new ResType(
					this,
					this.resTypes.size(),
					resType
				));
			}

			// index them by name
			resTypesByName = new HashMap<>();
			for (ResType rt : this.resTypes) {
				resTypesByName.put(rt.name, rt);
			}

			// get the wild type, if present for this position
			this.wildType = getResType(wildType);

			// gather the mutants
			mutations = this.resTypes.stream()
				.filter(rt -> rt != this.wildType)
				.collect(Collectors.toList());
		}

		public ResType getResType(String name) {
			return resTypesByName.get(name);
		}

		public ResType getResTypeOrThrow(String name) {
			ResType rt = getResType(name);
			if (rt != null) {
				return rt;
			}
			throw new NoSuchElementException("Res type " + name + " not allowed at position " + resNum + ". Try one of " + resTypes);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Position && equals((Position)other);
		}

		public boolean equals(Position other) {
			return this.index == other.index
				&& this.resNum.equals(other.resNum)
				&& this.wildType.equals(other.wildType)
				&& this.resTypes.equals(other.resTypes);
		}

		@Override
		public String toString() {
			return String.format("%d:%s", index, resNum);
		}

		@Override
		public int compareTo(Position other) {
			return this.index - other.index;
		}
	}

	public class ResType implements Comparable<ResType>, Serializable {

		private static final long serialVersionUID = 1628630395970846143L;

		public final Position pos;
		public final int index;
		public final String name;

		public ResType(Position pos, int index, String name) {
			this.pos = pos;
			this.index = index;
			this.name = name;
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof ResType && equals((ResType)other);
		}

		public boolean equals(ResType other) {
			return this.index == other.index
				&& this.name.equals(other.name);
		}

		@Override
		public String toString() {
			return String.format("%d:%s", index, name);
		}

		@Override
		public int compareTo(ResType other) {
			return this.index - other.index;
		}

		public boolean isWildType() {
			return pos.wildType == this;
		}

		public boolean isMutation() {
			return pos.wildType != this;
		}

		/**
		 * Returns the residue type name in lower case if wild-type, upper case if mutation.
		 */
		public String mutationName() {
			if (isWildType()) {
				return name.toLowerCase();
			} else {
				return name.toUpperCase();
			}
		}
	}

	public final List<Position> positions = new ArrayList<>();

	private final Map<String,Position> positionsByResNum = new HashMap<>();

	public SeqSpace(SimpleConfSpace confSpace) {
		for (SimpleConfSpace.Position pos : confSpace.mutablePositions) {
			makePos(pos.resNum, pos.resFlex.wildType, pos.resTypes);
		}
	}

	private SeqSpace() {
		// used by union()
	}

	private void makePos(String resNum, String wildType, List<String> resTypes) {
		Position pos = new Position(positions.size(), resNum, wildType, resTypes);
		positions.add(pos);
		positionsByResNum.put(resNum, pos);
	}

	public List<String> getResNums() {
		return positions.stream()
			.map(pos -> pos.resNum)
			.collect(Collectors.toList());
	}

	public Position getPosition(String resNum) {
		return positionsByResNum.get(resNum);
	}

	public Position getPositionOrThrow(String resNum) {
		Position pos = getPosition(resNum);
		if (pos != null) {
			return pos;
		}
		throw new NoSuchElementException("no pos found in this sequence space at residue " + resNum + ". try one of " + getResNums());
	}

	public boolean containsWildTypeSequence() {
		for (Position pos : positions) {
			if (pos.wildType == null) {
				return false;
			}
		}
		return true;
	}

	public Sequence makeUnassignedSequence() {
		return new Sequence(this);
	}

	public Sequence makeWildTypeSequence() {
		if (!containsWildTypeSequence()) {
			throw new NoSuchElementException("sequence space does not contain the wild-type sequence, so cannot create it");
		}
		Sequence seq = new Sequence(this);
		seq.fillWildType();
		return seq;
	}

	/**
	 * Make a sequence with the given residue type names
	 *
	 * residue type names must be given in the same order as the positions in the sequence space
	 */
	public Sequence makeSequence(List<String> resTypes) {

		// just in case...
		if (resTypes.size() != positions.size()) {
			throw new IllegalArgumentException(String.format("expected %d residue types, but only got %d: %s",
				positions.size(),
				resTypes.size(),
				resTypes
			));
		}

		Sequence seq = makeUnassignedSequence();
		for (int i=0; i<positions.size(); i++) {
			Position pos = positions.get(i);
			seq.set(pos, resTypes.get(i));
		}
		return seq;
	}

	public List<Sequence> getMutants() {
		return getMutants(positions.size());
	}

	public List<Sequence> getMutants(int maxSimultaneousMutations) {
		return getMutants(maxSimultaneousMutations, false);
	}

	public List<Sequence> getMutants(int maxSimultaneousMutations, boolean reversePositionOrder) {

		List<Sequence> sequences = new ArrayList<>();

		// get all possible combinations of mutations
		List<List<SeqSpace.Position>> powersetOfPositions = MathTools.powersetUpTo(positions, maxSimultaneousMutations);

		// reverse to match order of old K* code if needed
		if (reversePositionOrder) {
			Collections.reverse(powersetOfPositions);
		}

		for (List<SeqSpace.Position> positions : powersetOfPositions) {

			// collect the mutations (res types except for wild type) for these positions into a simple list list
			List<List<ResType>> mutationsByPos = positions.stream()
				.map(pos -> pos.mutations)
				.collect(Collectors.toList());

			// enumerate all the combinations of res types
			for (List<ResType> mutations : MathTools.cartesianProduct(mutationsByPos)) {

				// build the sequence and add the mutations
				Sequence sequence = makeUnassignedSequence();
				sequence.fillWildType();
				for (ResType rt : mutations) {
					sequence.set(rt.pos, rt);
				}

				// if we have a full sequence, output it
				if (sequence.isFullyAssigned()) {
					sequences.add(sequence);
				}

				// NOTE: if sequence is not fully assigned, that means we're missing wild types at some positions,
				// and we hit our mutation limits, so don't output this sequence
			}
		}

		return sequences;
	}

	@Override
	public String toString() {
		StringBuilder buf = new StringBuilder();
		buf.append("Residue Types:");
		for (SeqSpace.Position pos : positions) {
			buf.append("\n\t");
			buf.append(pos.index);
			buf.append(":");
			buf.append(pos.resNum);
			buf.append("  [");
			for (ResType rt : pos.resTypes) {
				buf.append(" ");
				buf.append(rt.index);
				buf.append(":");
				if (pos.wildType == rt) {
					buf.append(rt.name.toLowerCase());
				} else {
					buf.append(rt.name);
				}
			}
			buf.append(" ]");
		}
		return buf.toString();
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof SeqSpace && equals((SeqSpace)other);
	}

	public boolean equals(SeqSpace other) {
		return this.positions.equals(other.positions);
	}
}
