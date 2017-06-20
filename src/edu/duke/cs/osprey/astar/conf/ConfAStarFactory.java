package edu.duke.cs.osprey.astar.conf;

import edu.duke.cs.osprey.externalMemory.Queue;

public interface ConfAStarFactory {
	
	Queue<ConfAStarNode> makeQueue(RCs rcs);
	ConfAStarNode makeRootNode(int numPos);
}
