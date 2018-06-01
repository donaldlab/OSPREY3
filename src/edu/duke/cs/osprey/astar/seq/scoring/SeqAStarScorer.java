package edu.duke.cs.osprey.astar.seq.scoring;


import edu.duke.cs.osprey.astar.seq.SeqAStarNode;

public interface SeqAStarScorer {

	/**
	 * calculate a score for the given node
	 */
	double calc(SeqAStarNode node);

	/**
	 * calculate a score for a given assignment, using the parent node to optimize, if possible
	 */
	default double calcDifferential(SeqAStarNode node, int nextPos, int nextRt) {
		// by default, make a new node and use the regular scoring
		return calc(node.assign(nextPos, nextRt));
	}
}
