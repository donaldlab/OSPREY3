package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.List;

public class SequentialOrderer implements ConfRanker.Orderer {

	@Override
	public int getNextPosition(ConfRanker ranker, ConfIndex confIndex, RCs rcs, double queryScore) {

		// always just get the next unassigned position
		return confIndex.undefinedPos[0];
	}
}
