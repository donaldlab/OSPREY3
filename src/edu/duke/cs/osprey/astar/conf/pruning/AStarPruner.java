package edu.duke.cs.osprey.astar.conf.pruning;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;

public interface AStarPruner {

	/**
	 * check to see if a node should be pruned, after it comes off the heap
	 */
	boolean isPruned(ConfAStarNode node);

	/**
	 * check to see if a child node should be pruned before it's created and scored
	 */
	boolean isPruned(ConfAStarNode node, int nextPos, int nextRc);
}
