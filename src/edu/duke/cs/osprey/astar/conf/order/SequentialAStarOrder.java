package edu.duke.cs.osprey.astar.conf.order;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public class SequentialAStarOrder implements AStarOrder {
	
	// NOTE: there's basically no reason to ever use this except to benchmark better methods
	// use a more intelligent static ordering instead

	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		// don't care...
	}
	
	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		
		// easy peasy
		// eg, the root node has level 0, so expand pos 0 next
		return confIndex.node.getLevel();
	}
}
