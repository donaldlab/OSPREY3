package edu.duke.cs.osprey.astar.seq;

import edu.duke.cs.osprey.astar.OptimizableAStarNode;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;


public interface SeqAStarNode extends OptimizableAStarNode, Comparable<SeqAStarNode> {

	SeqAStarNode assign(int pos, int rt);
	void getSequence(Sequence seq);

	Object getData();
	void setData(Object data);

	@Override
	default int compareTo(SeqAStarNode other) {
		return Double.compare(this.getScore(), other.getScore());
	}

	default Sequence makeSequence(SimpleConfSpace confSpace) {
		Sequence seq = confSpace.makeUnassignedSequence();
		getSequence(seq);
		return seq;
	}
}
