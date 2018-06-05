package edu.duke.cs.osprey.astar.seq.scoring;


import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

public interface SeqAStarScorer {

	/**
	 * calculate a score for the given node
	 */
	double calc(SeqAStarNode.Assignments assignments);

	/**
	 * calculate a score for a given assignment, using the parent node to optimize, if possible
	 */
	default double calcDifferential(SeqAStarNode.Assignments assignments, int nextPos, int nextRt) {

		// by default update the assignments and punt to calc()
		assignments.assign(nextPos, nextRt);
		double score = calc(assignments);
		assignments.unassign(nextPos);
		return score;
	}
}
