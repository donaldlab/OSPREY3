package edu.duke.cs.osprey.astar.seq.scoring;

import edu.duke.cs.osprey.astar.seq.SeqAStarNode;

/**
 * ie, a no-op
 */
public class NOPSeqAStarScorer implements SeqAStarScorer {

	@Override
	public double calc(SeqAStarNode node) {
		return 0;
	}

	@Override
	public double calcDifferential(SeqAStarNode node, int nextPos, int nextRt) {
		return 0;
	}
}
