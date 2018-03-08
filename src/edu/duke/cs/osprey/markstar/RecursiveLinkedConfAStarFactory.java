package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.externalMemory.Queue;

public class RecursiveLinkedConfAStarFactory implements ConfAStarFactory {
	
	@Override
	public Queue<ConfAStarNode> makeQueue(RCs rcs) {
		
		Queue<RecursiveLinkedAStarNode> pq = Queue.PriorityFactory.of(null);
		
		// java's type system is dumb sometimes...
		Queue<? extends ConfAStarNode> q2 = (Queue<? extends ConfAStarNode>)pq;
		@SuppressWarnings("unchecked")
		Queue<ConfAStarNode> q3 = (Queue<ConfAStarNode>)q2;
		return q3;
	}
	
	@Override
	public ConfAStarNode makeRootNode(int numPos) {
		return new RecursiveLinkedAStarNode();
	}
	
	
	public RecursiveLinkedAStarNode makeRecRootNode(int numPos) {
		return new RecursiveLinkedAStarNode();
	}
	
	public Queue<RecursiveLinkedAStarNode> makeRecQueue(RCs rcs) {
		return Queue.PriorityFactory.of(null);
	}
}
