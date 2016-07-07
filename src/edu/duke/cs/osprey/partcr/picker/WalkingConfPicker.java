package edu.duke.cs.osprey.partcr.picker;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;

public class WalkingConfPicker implements ConfPicker {

	private int numItersPerConf;
	private Map<ConfAStarNode,Integer> iteratedNodes;
	
	public WalkingConfPicker() {
		this(1);
	}
	
	public WalkingConfPicker(int numItersPerConf) {
		this.numItersPerConf = numItersPerConf;
		iteratedNodes = new HashMap<>(); // yeah, ok to match on instance
	}

	@Override
	public ConfAStarNode pick(List<ConfAStarNode> nodes) {
		
		for (ConfAStarNode node : nodes) {
		
			// have we iterated this node/conf yet?
			Integer numIters = iteratedNodes.get(node);
			if (numIters == null) {
				
				// nope, let's start
				iteratedNodes.put(node, 1);
				return node;
				
			} else if (numIters < numItersPerConf) {
				
				// yup, but not enough
				iteratedNodes.put(node, numIters + 1);
				return node;
			}
			
			// yeah, but took much, so look to the next node
		}
		
		throw new IllegalStateException("ran out of conformations, can't pick one");
	}
}
