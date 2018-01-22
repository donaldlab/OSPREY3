package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.List;

public class SequentialOrderer implements ConfRanker.Orderer {

	@Override
	public SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore) {

		// always just get the next position
		return unassignedPositions.get(0);
	}
}
