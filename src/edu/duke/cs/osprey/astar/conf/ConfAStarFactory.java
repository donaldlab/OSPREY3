package edu.duke.cs.osprey.astar.conf;

import edu.duke.cs.osprey.externalMemory.Queue;

public interface ConfAStarFactory {
	
	Queue<ConfAStarNode> makeQueue();
	ConfAStarNode makeRootNode(int numPos);
}
