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

package edu.duke.cs.osprey.astar.seq.nodes;

import edu.duke.cs.osprey.astar.OptimizableAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.Arrays;


public interface SeqAStarNode extends OptimizableAStarNode, Comparable<SeqAStarNode> {

	SeqAStarNode assign(int pos, int rt);
	void getAssignments(Assignments assignments);
	void getSequence(Sequence seq);

	Object getData();
	void setData(Object data);

	@Override
	default int compareTo(SeqAStarNode other) {
		return Double.compare(this.getScore(), other.getScore());
	}

	default Sequence makeSequence(SeqSpace seqSpace) {
		Sequence seq = seqSpace.makeUnassignedSequence();
		getSequence(seq);
		return seq;
	}


	/**
	 * efficiently-accessible assignment info regardless of how nodes are implemented internally
	 *
	 * position orders must always be sorted in increasing order
	 */
	public static class Assignments {

		public final int numPos;

		public final int[] assignedPos;
		public final int[] assignedRTs;
		public int numAssigned;
		public final int[] unassignedPos;
		public int numUnassigned;

		public Assignments(int numPos) {
			this.numPos = numPos;
			assignedPos = new int[numPos];
			assignedRTs = new int[numPos];
			numAssigned = 0;
			unassignedPos = new int[numPos];
			for (int i = 0; i< numPos; i++) {
				unassignedPos[i] = i;
			}
			numUnassigned = numPos;
		}

		/**
		 * ensures assigned and unassigned positions are sorted in increasing order
		 */
		public void sortAssigned() {

			// these arrays are typically tiny (<20), so insertion sort is very efficient here
			// also, we have to sort the assigned positions and res types together,
			// so we can't use a library sort =(

			// sort the asssigned side
			for (int i=1; i<numAssigned; i++) {

				int tempPos = assignedPos[i];
				int tempRT = assignedRTs[i];

				int j;
				for (j=i; j>=1 && tempPos < assignedPos[j-1]; j--) {
					assignedPos[j] = assignedPos[j-1];
					assignedRTs[j] = assignedRTs[j-1];
				}
				assignedPos[j] = tempPos;
				assignedRTs[j] = tempRT;
			}
		}

		public Integer getAssignment(int pos) {
			int i = Arrays.binarySearch(assignedPos, 0, numAssigned, pos);
			if (i >= 0) {
				return assignedRTs[i];
			}
			return null;
		}

		public boolean isAssigned(int pos) {
			return getAssignment(pos) != null;
		}

		/**
		 * Populates the unassigned positions, based on what's not assigned
		 * @return this, for method chaining
		 */
		public Assignments updateUnassigned() {
			numUnassigned = 0;
			for (int pos=0; pos<numPos; pos++) {
				if (!isAssigned(pos)) {
					unassignedPos[numUnassigned++] = pos;
				}
			}
			return this;
		}

		/**
		 * Assigns the residue type to the position
		 * @return this, for method chaining
		 */
		public Assignments assign(int pos, int rt) {

			// find the pos in the unassigned section
			int i = Arrays.binarySearch(unassignedPos, 0, numUnassigned, pos);
			if (i < 0) {
				throw new IllegalArgumentException("position " + pos + " is already assigned");
			}

			// remove the unassigned pos
			numUnassigned--;
			System.arraycopy(unassignedPos, i + 1, unassignedPos, i, numUnassigned - i);

			// make room for the new assignment
			i = Arrays.binarySearch(assignedPos, 0, numAssigned, pos);
			assert (i < 0); // the pos shouldn't not be assigned
			i = -i - 1; // transform into insertion pos
			System.arraycopy(assignedPos, i, assignedPos, i + 1, numAssigned - i);
			System.arraycopy(assignedRTs, i, assignedRTs, i + 1, numAssigned - i);
			numAssigned++;

			// finally, make the assignment
			assignedPos[i] = pos;
			assignedRTs[i] = rt;

			return this;
		}

		/**
		 * Removes the assignment from the position
		 * @return this, for method chaining
		 */
		public Assignments unassign(int pos) {

			// find the pos in the assigned section
			int i = Arrays.binarySearch(assignedPos, 0, numAssigned, pos);
			if (i < 0) {
				throw new IllegalArgumentException("position " + pos + " is not assigned");
			}

			// remove the assigned pos
			numAssigned--;
			System.arraycopy(assignedPos, i + 1, assignedPos, i, numAssigned - i);
			System.arraycopy(assignedRTs, i + 1, assignedRTs, i, numAssigned - i);

			// make room for the new unassignment
			i = Arrays.binarySearch(unassignedPos, 0, numUnassigned, pos);
			assert (i < 0); // the pos shouldn't not be assigned
			i = -i - 1; // transform into insertion pos
			System.arraycopy(unassignedPos, i, unassignedPos, i + 1, numUnassigned - i);
			numUnassigned++;

			unassignedPos[i] = pos;

			return this;
		}

		public RCs makeRCs(SeqSpace seqSpace, SimpleConfSpace confSpace) {
			return new RCs(confSpace, (confPos, resConf) -> {

				// immutable pos? keep everything
				if (!confPos.hasMutations()) {
					return true;
				}

				// map the conf space pos to the seq space pos
				SeqSpace.Position seqPos = seqSpace.getPosition(confPos.resNum);

				// is this seq pos unassigned? keep everything
				Integer assignment = getAssignment(seqPos.index);
				if (assignment == null) {
					return true;
				}

				// otherwise, only keep matching RCs
				SeqSpace.ResType resType = seqPos.resTypes.get(assignment);
				return resType.name.equals(resConf.template.name);
			});
		}

		@Override
		public String toString() {
			StringBuilder buf = new StringBuilder();
			for (int i=0; i<numAssigned; i++) {
				if (buf.length() > 0) {
					buf.append(" ");
				}
				buf.append(assignedPos[i]);
				buf.append("=");
				buf.append(assignedRTs[i]);
			}
			return buf.toString();
		}

		public String toString(SeqSpace seqSpace) {
			StringBuilder buf = new StringBuilder();
			for (int i=0; i<numAssigned; i++) {
				SeqSpace.Position pos = seqSpace.positions.get(assignedPos[i]);
				SeqSpace.ResType rt = pos.resTypes.get(assignedRTs[i]);
				if (buf.length() > 0) {
					buf.append(" ");
				}
				buf.append(pos.toString());
				buf.append("=");
				buf.append(rt.toString());
			}
			return buf.toString();
		}
	}
}
