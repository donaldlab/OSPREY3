package edu.duke.cs.osprey.astar.seq.scoring;

import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

/**
 * ie, a no-op
 */
public class NOPSeqAStarScorer implements SeqAStarScorer {

	@Override
	public double calc(SeqAStarNode.Assignments assignments) {
		return 0;
	}

	@Override
	public double calcDifferential(SeqAStarNode.Assignments assignments, int nextPos, int nextRt) {
		return 0;
	}
}
