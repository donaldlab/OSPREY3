package edu.duke.cs.osprey.astar.conf.order;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public class SequentialAStarOrder implements AStarOrder {

	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		// don't care...
	}
	
	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		
		// easy peasy
		// eg, the root node has level 0, so expand pos 0 next
		return confIndex.getNode().getLevel();
	}
}
