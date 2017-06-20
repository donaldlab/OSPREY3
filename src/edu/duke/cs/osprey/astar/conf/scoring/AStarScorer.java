package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;

public interface AStarScorer {

	AStarScorer make();
	double calc(ConfIndex confIndex, RCs rcs);
	
	default double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
		// just punt to calc() by default
		return calc(confIndex.assign(nextPos, nextRc), rcs);
	}
}
