package edu.duke.cs.osprey.astar.conf.order;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public interface AStarOrder {

	void setScorers(AStarScorer gscorer, AStarScorer hscorer);
	int getNextPos(ConfIndex confIndex, RCs rcs);
}
