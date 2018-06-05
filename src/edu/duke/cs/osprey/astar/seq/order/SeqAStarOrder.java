package edu.duke.cs.osprey.astar.seq.order;

import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;

public interface SeqAStarOrder {
	int getNextPos(SeqAStarNode node, RTs rts);
}
