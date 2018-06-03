package edu.duke.cs.osprey.astar.seq;

import edu.duke.cs.osprey.astar.OptimizableAStarNode;
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

	default Sequence makeSequence(SimpleConfSpace confSpace) {
		Sequence seq = confSpace.makeUnassignedSequence();
		getSequence(seq);
		return seq;
	}


	/**
	 * efficiently-accessible assignment info regardless of how nodes are implemented internally
	 *
	 * position orders must always be sorted in increasing order
	 */
	public static class Assignments {

		public final int numMutablePos;

		public final int[] assignedMPos;
		public final int[] assignedRTs;
		public int numAssigned;
		public final int[] unassignedMPos;
		public int numUnassigned;

		public Assignments(int numMutablePos) {
			this.numMutablePos = numMutablePos;
			assignedMPos = new int[numMutablePos];
			assignedRTs = new int[numMutablePos];
			numAssigned = 0;
			unassignedMPos = new int[numMutablePos];
			for (int i=0; i<numMutablePos; i++) {
				unassignedMPos[i] = i;
			}
			numUnassigned = numMutablePos;
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

				int tempPos = assignedMPos[i];
				int tempRT = assignedRTs[i];

				int j;
				for (j=i; j>=1 && tempPos < assignedMPos[j-1]; j--) {
					assignedMPos[j] = assignedMPos[j-1];
					assignedRTs[j] = assignedRTs[j-1];
				}
				assignedMPos[j] = tempPos;
				assignedRTs[j] = tempRT;
			}

			// sort the unassigned side
			for (int i=1; i<numUnassigned; i++) {

				int tempPos = unassignedMPos[i];

				int j;
				for (j=i; j>=1 && tempPos < unassignedMPos[j-1]; j--) {
					unassignedMPos[j] = unassignedMPos[j-1];
				}
				unassignedMPos[j] = tempPos;
			}
		}

		public Integer getAssignment(int mpos) {
			int i = Arrays.binarySearch(assignedMPos, 0, numAssigned, mpos);
			if (i >= 0) {
				return assignedRTs[i];
			}
			return null;
		}

		public boolean isAssigned(int mpos) {
			return getAssignment(mpos) != null;
		}
	}
}
