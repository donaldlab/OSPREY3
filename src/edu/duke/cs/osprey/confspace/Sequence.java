package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.tools.HashCalculator;
import org.apache.commons.lang3.StringUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

public class Sequence implements Iterable<Sequence.Assignment> {

	public class Assignment {

		private SimpleConfSpace.Position pos;

		public Assignment(SimpleConfSpace.Position pos) {
			this.pos = pos;
		}

		/**
		 * Get the position in the conformation space of this assignment.
		 */
		public SimpleConfSpace.Position getPos() {
			return pos;
		}

		/**
		 * Get the residue number of this assignment.
		 */
		public String getResNum() {
			return pos.resNum;
		}

		/**
		 * Get the residue type of this assignment.
		 */
		public String getResType() {
			return resTypes[pos.index];
		}

		/**
		 * Set the residue type of this assignment.
		 */
		public void setResType(String resType) {
			set(pos, resType);
		}

		/**
		 * Get the residue type of the wild-type at the residue of this assignment.
		 */
		public String getWildType() {
			return Sequence.this.getWildType(pos);
		}

		/**
		 * Set the residue type of this assignment to the wild-type.
		 */
		public void setWildType() {
			Sequence.this.setWildType(pos);
		}

		/**
		 * Return true if this residue is assigned
		 */
		public boolean isAssigned() {
			return Sequence.this.isAssigned(pos);
		}

		/**
		 * Return true if this residue is assigned and matches the wild-type.
		 */
		public boolean isWildType() {
			return Sequence.this.isWildType(pos);
		}

		/**
		 * Return true if this residue is assigned and differs from the wild-type.
		 */
		public boolean isMutated() {
			return Sequence.this.isMutated(pos);
		}

		@Override
		public String toString() {
			return String.format("%s=%s", getResNum(), getResType());
		}
	}

	/**
	 * The conformation space associated with this sequence
	 */
	public final SimpleConfSpace confSpace;

	private final String[] resTypes;

	/**
	 * Make a sequence with no assignments.
	 */
	public static Sequence makeUnassigned(SimpleConfSpace confSpace) {
		return new Sequence(confSpace);
	}

	/**
	 * Make a sequence set to the wild-type residue types.
	 */
	public static Sequence makeWildType(SimpleConfSpace confSpace) {
		return makeUnassigned(confSpace).fillWildType();
	}

	/**
	 * Make a sequence from assignments to design positions
	 */
	public static Sequence makeFromAssignments(SimpleConfSpace confSpace, int[] assignments) {
		Sequence sequence = makeUnassigned(confSpace);
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			sequence.set(pos, pos.resConfs.get(assignments[pos.index]).template.name);
		}
		return sequence;
	}

	/**
	 * Make a sequence from a conformation
	 */
	public static Sequence makeFromConf(SimpleConfSpace confSpace, ConfSearch.ScoredConf conf) {
		return makeFromAssignments(confSpace, conf.getAssignments());
	}

	private Sequence(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
		this.resTypes = new String[confSpace.positions.size()];
		Arrays.fill(resTypes, null);
	}

	/**
	 * Copy constructor.
	 */
	public Sequence(Sequence other) {
		this.confSpace = other.confSpace;
		this.resTypes = other.resTypes.clone();
	}

	/**
	 * Get the residue type at the specified position.
	 */
	public String get(SimpleConfSpace.Position pos) {
		return resTypes[pos.index];
	}

	/**
	 * Get the residue type at the specified residue.
	 */
	public String get(String resNum) {
		return get(confSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Set the residue type at the specified position.
	 * @return this sequence, for method chaining
	 */
	public Sequence set(SimpleConfSpace.Position pos, String resType) {
		resTypes[pos.index] = resType;
		return this;
	}

	/**
	 * Set the residue type at the specified residue.
	 * @return this sequence, for method chaining
	 */
	public Sequence set(String resNum, String resType) {
		set(confSpace.getPositionOrThrow(resNum), resType);
		return this;
	}

	private class AssignmentIterator implements Iterator<Assignment> {

		private Iterator<SimpleConfSpace.Position> posIter;

		public AssignmentIterator(Iterable<SimpleConfSpace.Position> positions) {
			posIter = positions.iterator();
		}

		@Override
		public boolean hasNext() {
			return posIter.hasNext();
		}

		@Override
		public Assignment next() {
			return new Assignment(posIter.next());
		}
	}

	public Iterator<Assignment> iterator() {
		return iterator(confSpace.positions);
	}

	public Iterator<Assignment> iterator(Iterable<SimpleConfSpace.Position> positions) {
		return new AssignmentIterator(positions);
	}

	public Stream<Assignment> stream() {
		return stream(confSpace.positions);
	}

	public Stream<Assignment> stream(Iterable<SimpleConfSpace.Position> positions) {
		return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(positions), 0), false);
	}

	/**
	 * Returns the residue type of the wild-type at the specified position
	 */
	public String getWildType(SimpleConfSpace.Position pos) {
		return pos.resFlex.wildType;
	}

	/**
	 * Returns the residue type of the wild-type at the specified residue
	 */
	public String getWildType(String resNum) {
		return getWildType(confSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Sets the residue type at the specified position to wild-type
	 * @return this sequence, for method chaining
	 */
	public Sequence setWildType(SimpleConfSpace.Position pos) {
		set(pos, getWildType(pos));
		return this;
	}

	/**
	 * Sets the residue type at the specified residue to wild-type
	 * @return this sequence, for method chaining
	 */
	public Sequence setWildType(String resNum) {
		setWildType(confSpace.getPositionOrThrow(resNum));
		return this;
	}

	/**
	 * Set all unassigned residues to the wild-type residue type.
	 * @return this sequence, for method chaining
	 */
	public Sequence fillWildType() {
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (get(pos) == null) {
				set(pos, getWildType(pos));
			}
		}
		return this;
	}

	/**
	 * Make a new sequence with the mutation at the specified position.
	 */
	public Sequence makeMutatedSequence(SimpleConfSpace.Position pos, String resType) {
		return new Sequence(this)
			.set(pos, resType);
	}

	/**
	 * Make a new sequence with the mutation at the specified residue.
	 */
	public Sequence makeMutatedSequence(String resNum, String resType) {
		return makeMutatedSequence(confSpace.getPositionOrThrow(resNum), resType);
	}

	/**
	 * Return true if the specified position is assigned.
	 */
	public boolean isAssigned(SimpleConfSpace.Position pos) {
		return get(pos) != null;
	}

	/**
	 * Return true if the specified position is assigned.
	 */
	public boolean isAssigned(String resNum) {
		return isAssigned(confSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Count the number of assigned residues.
	 */
	public int countAssignments() {
		int count = 0;
		for (String resType : resTypes) {
			if (resType != null) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Return true if all residues are assigned, false if not.
	 */
	public boolean isFullyAssigned() {
		for (String resType : resTypes) {
			if (resType == null) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true if the specified position is assigned, and the residue type matches the wild-type.
	 */
	public boolean isWildType(SimpleConfSpace.Position pos) {
		return isAssigned(pos) && get(pos).equalsIgnoreCase(pos.resFlex.wildType);
	}

	/**
	 * Return true if the specified residue is assigned, and the residue type matches the wild-type.
	 */
	public boolean isWildType(String resNum) {
		return isWildType(confSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Return true if the residue types of all residues match the wild-type.
	 */
	public boolean isWildType() {
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (!isWildType(pos)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Return true if the specified position is assigned, and the residue type differs from the wild-type.
	 */
	public boolean isMutated(SimpleConfSpace.Position pos) {
		return isAssigned(pos) && !get(pos).equalsIgnoreCase(pos.resFlex.wildType);
	}

	/**
	 * Return true if the specified residue is assigned, and the residue type differs from the wild-type.
	 */
	public boolean isMutated(String resNum) {
		return isMutated(confSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Count the number of residues whose residue types are different from the wild-type.
	 * Ignores unassigned residues.
	 */
	public int countMutations() {
		int count = 0;
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (isMutated(pos)) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Make a new sequence that is a subset of this sequence, using the positions in the specified conformation space.
	 * The positions in the specified conformation space must exist in this sequence's conformation space.
	 */
	public Sequence filter(SimpleConfSpace confSpace) {
		Sequence subsequence = confSpace.makeUnassignedSequence();
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			subsequence.set(pos.resNum, this.get(pos.resNum));
		}
		return subsequence;
	}

	public RCs makeRCs() {
		return new RCs(confSpace, (pos, resConf) -> {

			// no assignments? keep everything
			if (!isAssigned(pos)) {
				return true;
			}

			// otherwise, only keep matching RCs
			return get(pos).equals(resConf.template.name);
		});
	}

	@Override
	public String toString() {
		return toString(Renderer.Assignment, null, confSpace.positions);
	}

	public static enum Renderer {

		/**
		 * Render just the residue number.
		 */
		ResNum {
			@Override
			public String render(Sequence.Assignment assignment) {
				return assignment.getResNum();
			}
		},

		/**
		 * Render just the residue type.
		 */
		ResType {
			@Override
			public String render(Sequence.Assignment assignment) {
				return assignment.getResType();
			}
		},

		/**
		 * Render just the residue type, but use UPPERCASE for mutations, and lowercase for wild-type.
		 */
		ResTypeMutations {
			@Override
			public String render(Sequence.Assignment assignment) {
				if (assignment.isMutated()) {
					return assignment.getResType().toUpperCase();
				} else {
					return assignment.getResType().toLowerCase();
				}
			}
		},

		/**
		 * Render "R=T" where R is the residue number, and T is the residue type.
		 */
		Assignment {
			@Override
			public String render(Sequence.Assignment assignment) {
				return ResNum.render(assignment) + "=" + ResType.render(assignment);
			}
		},

		/**
		 * Render "R=T" where R is the residue number, and T is the residue type,
		 * but use UPPERCASE for mutations, and lowercase for wild-type.
		 */
		AssignmentMutations {
			@Override
			public String render(Sequence.Assignment assignment) {
				return ResNum.render(assignment) + "=" + ResTypeMutations.render(assignment);
			}
		};

		public abstract String render(Assignment assignment);
	}

	public String toString(Renderer renderer) {
		return toString(renderer, null, confSpace.positions);
	}

	public String toString(Renderer renderer, int cellSize) {
		return toString(renderer, cellSize, confSpace.positions);
	}

	/**
	 * Render the sequence into a string.
	 * @param renderer How to render each assignment.
	 * @param cellSize Right pad each assignment to this amount, or null for no padding. Useful for horizontal alignment.
	 * @param positions Render assignments in order of these positions.
	 * @return A string representing the sequence.
	 */
	public String toString(Renderer renderer, Integer cellSize, List<SimpleConfSpace.Position> positions) {
		return String.join(" ", stream(positions)
			.map((assignment) -> {

				// run the assignment renderer
				String rendered = renderer.render(assignment);

				// apply cell size if needed
				if (cellSize != null) {
					rendered = StringUtils.rightPad(rendered, cellSize, ' ');
				}

				return rendered;
			})
			.collect(Collectors.toList())
		);
	}

	public int getMaxResNumLength() {
		return confSpace.positions.stream()
			.map((pos) -> pos.resNum.length())
			.max(Integer::compare)
			.get();
	}

	private String[] getResNums() {
		return confSpace.positions.stream()
			.map((pos) -> pos.resNum)
			// apparently python doc tool 'javalang' can't parse this
			//.toArray(String[]::new);
			// so use old-fashioned way
			.toArray((arraySize) -> new String[arraySize]);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof Sequence && equals((Sequence)other);
	}

	public boolean equals(Sequence other) {
		return Arrays.equals(this.resTypes, other.resTypes)
			&& Arrays.equals(this.getResNums(), other.getResNums());
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineHashes(
			Arrays.hashCode(resTypes),
			Arrays.hashCode(getResNums())
		);
	}
}
