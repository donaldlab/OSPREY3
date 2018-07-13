package edu.duke.cs.osprey.astar.seq.order;

import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

public class SequentialSeqAStarOrder implements SeqAStarOrder {
	
	@Override
	public int getNextPos(SeqAStarNode.Assignments assignments, RTs rts) {
		
		// easy peasy, pick the first unassigned pos
		if (assignments.numUnassigned > 0) {
			return assignments.unassignedPos[0];
		}

		throw new IllegalArgumentException("node has no unassigned positions");
	}
}
