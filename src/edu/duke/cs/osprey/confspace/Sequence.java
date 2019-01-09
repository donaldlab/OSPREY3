/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.Streams;
import org.apache.commons.lang3.StringUtils;

import java.util.*;
import java.util.stream.Collectors;


public class Sequence {

	public class Assignment {

		public final SeqSpace.Position pos;

		private Assignment(SeqSpace.Position pos) {
			this.pos = pos;
		}

		/**
		 * Get the residue number to which this assignment refers
		 */
		public String getResNum() {
			return pos.resNum;
		}

		/**
		 * Get the residue type of this assignment, or null if none
		 */
		public SeqSpace.ResType getResType() {
			return Sequence.this.get(pos);
		}

		/**
		 * Set the residue type of this assignment.
		 */
		public void setResType(String resType) {
			Sequence.this.set(pos, resType);
		}

		/**
		 * Get the residue type of the wild-type at the residue of this assignment.
		 */
		public SeqSpace.ResType getWildType() {
			return pos.wildType;
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
	}

	public static final int Unassigned = -1;

	public final SeqSpace seqSpace;

	public final int[] rtIndices;

	public Sequence(SeqSpace seqSpace) {
		this.seqSpace = seqSpace;
		rtIndices = new int[seqSpace.positions.size()];
		Arrays.fill(rtIndices, Unassigned);
	}

	public Sequence(SeqSpace seqSpace, int[] rtIndices) {
		this.seqSpace = seqSpace;
		this.rtIndices = rtIndices;
	}

	/**
	 * Returns the list of residue numbers for this sequence.
	 */
	public List<String> resNums() {
		return seqSpace.positions.stream()
			.map(pos -> pos.resNum)
			.collect(Collectors.toList());
	}

	private void checkPos(SeqSpace.Position pos) {
		if (pos.seqSpace != seqSpace) {
			throw new IllegalArgumentException("sequence position " + pos + " is not part of this sequence space");
		}
	}

	private void checkRT(SeqSpace.Position pos, SeqSpace.ResType rt) {
		if (rt.pos != pos) {
			throw new IllegalArgumentException("res type " + rt + " is not part of sequence position " + pos);
		}
	}

	/**
	 * Get the residue type at the specified residue by sequence position
	 */
	public SeqSpace.ResType get(SeqSpace.Position pos) {
		checkPos(pos);
		int rt = rtIndices[pos.index];
		if (rt == Unassigned) {
			return null;
		}
		return pos.resTypes.get(rt);
	}

	/**
	 * Get the residue type at the specified residue by residue number
	 */
	public SeqSpace.ResType get(String resNum) {
		return get(seqSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Set the residue type at the specified sequence position.
	 * @return this sequence, for method chaining
	 */
	public Sequence set(SeqSpace.Position pos, SeqSpace.ResType rt) {
		checkPos(pos);
		checkRT(pos, rt);
		rtIndices[pos.index] = rt.index;
		return this;
	}

	/**
	 * Set the residue type name at the specified residue number.
	 * @return this sequence, for method chaining
	 */
	public Sequence set(SeqSpace.Position pos, String resType) {
		checkPos(pos);
		SeqSpace.ResType rt = pos.getResTypeOrThrow(resType);
		rtIndices[pos.index] = rt.index;
		return this;
	}

	/**
	 * Set the residue type name at the specified residue number.
	 * @return this sequence, for method chaining
	 */
	public Sequence set(String resNum, String resType) {
		SeqSpace.Position pos = seqSpace.getPositionOrThrow(resNum);
		SeqSpace.ResType rt = pos.getResTypeOrThrow(resType);
		rtIndices[pos.index] = rt.index;
		return this;
	}

	/**
	 * Returns the residue type of the wild-type at the specified residue
	 */
	public SeqSpace.ResType getWildType(String resNum) {
		return seqSpace.getPositionOrThrow(resNum).wildType;
	}

	/**
	 * Sets the residue type at the specified sequence position to wild-type
	 * @return this sequence, for method chaining
	 */
	public Sequence setWildType(SeqSpace.Position pos) {
		checkPos(pos);
		if (pos.wildType == null) {
			throw new NoSuchElementException("sequence space does not contain the wild type at position " + pos);
		}
		rtIndices[pos.index] = pos.wildType.index;
		return this;
	}

	/**
	 * Sets the residue type at the specified residue to wild-type
	 * @return this sequence, for method chaining
	 */
	public Sequence setWildType(String resNum) {
		SeqSpace.Position pos = seqSpace.getPositionOrThrow(resNum);
		rtIndices[pos.index] = pos.wildType.index;
		return this;
	}

	/**
	 * Set all unassigned residues to the wild-type residue type, when the wild-type is part of the sequence space
	 * @return this sequence, for method chaining
	 */
	public Sequence fillWildType() {
		for (SeqSpace.Position pos : seqSpace.positions) {
			if (pos.wildType != null && !isAssigned(pos)) {
				setWildType(pos);
			}
		}
		return this;
	}

	public Iterable<Assignment> assignments() {
		return () -> seqSpace.positions.stream()
			.map(pos -> new Assignment(pos))
			.iterator();
	}

	public Iterable<Assignment> assignments(Iterable<String> resNums) {
		return () -> Streams.of(resNums)
			.map(resNum -> new Assignment(seqSpace.getPositionOrThrow(resNum)))
			.iterator();
	}

	/**
	 * Makes an identical copy of this sequence.
	 */
	public Sequence copy() {
		Sequence other = new Sequence(seqSpace);
		for (SeqSpace.Position pos : seqSpace.positions) {
			other.rtIndices[pos.index] = this.rtIndices[pos.index];
		}
		return other;
	}

	@Override
	public String toString() {
		return toString(Renderer.AssignmentMutations);
	}

	/**
	 * Return true if the specified sequence position is assigned.
	 */
	public boolean isAssigned(SeqSpace.Position pos) {
		return get(pos) != null;
	}

	/**
	 * Return true if the specified residue is assigned.
	 */
	public boolean isAssigned(String resNum) {
		return get(resNum) != null;
	}

	/**
	 * Count the number of assigned residues.
	 */
	public int countAssignments() {
		return (int)resNums().stream()
			.filter(resNum -> isAssigned(resNum))
			.count();
	}

	/**
	 * Return true if all residues are assigned, false if not.
	 */
	public boolean isFullyAssigned() {
		return resNums().stream()
			.allMatch(resNum -> isAssigned(resNum));
	}

	/**
	 * Return true if the specified sequence position is assigned, and the residue type matches the wild-type.
	 */
	public boolean isWildType(SeqSpace.Position pos) {
		checkPos(pos);
		SeqSpace.ResType resType = get(pos);
		return resType != null && resType == pos.wildType;
	}

	/**
	 * Return true if the specified residue is assigned, and the residue type matches the wild-type.
	 */
	public boolean isWildType(String resNum) {
		return isWildType(seqSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Return true if the residue types of all residues match the wild-type.
	 */
	public boolean isWildType() {
		return resNums().stream()
			.allMatch(resNum -> isWildType(resNum));
	}

	/**
	 * Return true if the specified seqeucence position is assigned, and the residue type differs from the wild-type.
	 */
	public boolean isMutated(SeqSpace.Position pos) {
		SeqSpace.ResType resType = get(pos);
		return resType != null && resType != pos.wildType;
	}

	/**
	 * Return true if the specified residue is assigned, and the residue type differs from the wild-type.
	 */
	public boolean isMutated(String resNum) {
		return isMutated(seqSpace.getPositionOrThrow(resNum));
	}

	/**
	 * Count the number of residues whose residue types are different from the wild-type.
	 * Ignores unassigned residues.
	 */
	public int countMutations() {
		return (int)resNums().stream()
			.filter(resNum -> isMutated(resNum))
			.count();
	}

	/**
	 * Returns a new sequence containing the subset of mutable positions belonging to the smaller sequence space.
	 */
	public Sequence filter(SeqSpace smallerSeqSpace) {
		Sequence seq = smallerSeqSpace.makeUnassignedSequence();
		for (SeqSpace.Position smallerPos : smallerSeqSpace.positions) {
			SeqSpace.Position biggerPos = seqSpace.getPositionOrThrow(smallerPos.resNum);
			SeqSpace.ResType biggerRT = get(biggerPos);
			if (biggerRT == null)
				System.out.println("PROBLEM");
			SeqSpace.ResType smallerRT = smallerPos.getResTypeOrThrow(biggerRT.name);
			seq.set(smallerPos, smallerRT);
		}
		return seq;
	}

	public RCs makeRCs(SimpleConfSpace confSpace) {
		return new RCs(confSpace, (pos, resConf) -> {

			// immutable pos? keep everything
			if (!pos.hasMutations()) {
				return true;
			}

			// no assignments? keep everything
			if (!isAssigned(pos.resNum)) {
				return true;
			}

			// otherwise, only keep matching RCs
			return get(pos.resNum).name.equals(resConf.template.name);
		});
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineHashes(
			seqSpace.hashCode(),
			Arrays.hashCode(rtIndices)
		);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof Sequence && equals((Sequence)other);
	}

	public boolean equals(Sequence other) {
		return this.seqSpace == other.seqSpace
			&& Arrays.equals(this.rtIndices, other.rtIndices);
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
				return assignment.getResType().name;
			}
		},

		/**
		 * Render just the residue type, but use UPPERCASE for mutations, and lowercase for wild-type.
		 */
		ResTypeMutations {
			@Override
			public String render(Sequence.Assignment assignment) {
				return assignment.getResType().mutationName();
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
		return toString(renderer, null, resNums());
	}

	public String toString(Renderer renderer, int cellSize) {
		return toString(renderer, cellSize, resNums());
	}

	/**
	 * Render the sequence into a string.
	 * @param renderer How to render each assignment.
	 * @param cellSize Right pad each assignment to this amount, or null for no padding. Useful for horizontal alignment.
	 * @param resNums Render assignments in order of these residues.
	 * @return A string representing the sequence.
	 */
	public String toString(Renderer renderer, Integer cellSize, List<String> resNums) {

		if (resNums.isEmpty()) {
			return "(no mutable residues)";
		}

		return String.join(" ", Streams.of(assignments(resNums))
			.filter(assignment -> assignment.isAssigned())
			.map(assignment -> {

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

	public int calcCellSize() {
		return seqSpace.positions.stream()
			.map(pos ->
				Math.max(
					pos.resNum.length(),
					pos.resTypes.get(rtIndices[pos.index]).name.length()
				)
			)
			.max(Integer::compare)
			.orElse(1);
	}
}
