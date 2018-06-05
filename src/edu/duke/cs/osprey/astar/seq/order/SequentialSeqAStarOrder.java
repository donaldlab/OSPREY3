package edu.duke.cs.osprey.astar.seq.order;

import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

public class SequentialSeqAStarOrder implements SeqAStarOrder {
	
	@Override
	public int getNextPos(SeqAStarNode node, RTs rts) {
		
		// easy peasy
		// eg, the root node has level 0, so expand pos 0 next
		return node.getLevel();
	}
}
