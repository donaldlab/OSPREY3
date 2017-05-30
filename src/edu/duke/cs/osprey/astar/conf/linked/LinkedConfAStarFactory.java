package edu.duke.cs.osprey.astar.conf.linked;

import edu.duke.cs.osprey.astar.conf.ConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.externalMemory.Queue;

public class LinkedConfAStarFactory implements ConfAStarFactory {
	
	@Override
	public Queue<ConfAStarNode> makeQueue(RCs rcs) {
		return Queue.PriorityFactory.of(null);
	}
	
	@Override
	public LinkedConfAStarNode makeRootNode(int numPos) {
		return new LinkedConfAStarNode();
	}
}
