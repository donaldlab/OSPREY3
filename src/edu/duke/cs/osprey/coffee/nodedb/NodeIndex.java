package edu.duke.cs.osprey.coffee.nodedb;

import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;

import java.util.Arrays;
import java.util.Deque;


/**
 * A specialization of FixedIndex for COFFEE nodes.
 */
public class NodeIndex {

	public static class Node implements FixedIndex.Indexable<BigExp> {

		public final int statei;
		public final int[] conf;
		public final BigExp zSumUpper;
		public final BigExp score;

		public Node(int statei, int[] conf, BigExp zSumUpper, BigExp score) {
			this.statei = statei;
			this.conf = conf;
			this.zSumUpper = zSumUpper;
			this.score = score;
		}

		/** re-score a node */
		public Node(Node other, BigExp score) {
			this(other.statei, other.conf, other.zSumUpper, score);
		}

		@Override
		public BigExp score() {
			return score;
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Node && equals((Node)other);
		}

		public boolean equals(Node other) {
			return this.statei == other.statei
				&& Arrays.equals(this.conf, other.conf)
				&& this.zSumUpper.equals(other.zSumUpper)
				&& this.score.equals(other.score);
		}

		@Override
		public String toString() {
			return String.format("Node[%d,%s,%s,%s]", statei, Arrays.toString(conf), zSumUpper, score);
		}

		public boolean isLeaf() {
			return Conf.isCompletelyAssigned(conf);
		}
	}

	public final BlockStore store;
	public final MultiStateConfSpace.State state;

	private final FixedIndex<BigExp,Node> index;

	public NodeIndex(BlockStore store, MultiStateConfSpace.State state) {

		this.store = store;
		this.state = state;

		index = new FixedIndex<>(store, Serializers.indexNode(state));
	}

	protected int nodesPerBlock() {
		return index.blockCapacity;
	}

	public long size() {
		return index.size();
	}

	public long freeSpace() {
		return index.freeSpace();
	}

	public boolean add(Node node) {
		assert (node.statei == state.index);
		return index.add(node);
	}

	public BigExp highestScore() {
		return index.highestScore();
	}

	public Node removeHighest() {
		return index.removeHighest();
	}

	public void freeUpSpace() {
		index.freeUpSpace();
	}

	public Deque<Node> dropped() {
		return index.dropped;
	}

	public void clear() {
		index.clear();
	}
}
