package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;

/**
 * ie, a no-op
 */
public class NOPAStarScorer implements AStarScorer {

	@Override
	public AStarScorer make() {
		return new NOPAStarScorer();
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		return 0;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		return 0;
	}
}
