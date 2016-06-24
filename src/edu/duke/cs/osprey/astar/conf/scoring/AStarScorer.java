package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;

public interface AStarScorer {

	double calc(ConfIndex confIndex, RCs rcs);
	double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc);
}
