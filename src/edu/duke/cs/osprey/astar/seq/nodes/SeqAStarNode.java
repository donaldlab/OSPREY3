package edu.duke.cs.osprey.astar.seq.nodes;

import edu.duke.cs.osprey.astar.OptimizableAStarNode;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;

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
		public void sort() {

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

			// sort the unassigned side
			for (int i=1; i<numUnassigned; i++) {

				int tempPos = unassignedPos[i];

				int j;
				for (j=i; j>=1 && tempPos < unassignedPos[j-1]; j--) {
					unassignedPos[j] = unassignedPos[j-1];
				}
				unassignedPos[j] = tempPos;
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
	}
}
