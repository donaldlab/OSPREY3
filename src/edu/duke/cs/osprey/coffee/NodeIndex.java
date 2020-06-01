package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;

import java.util.Arrays;
import java.util.function.Consumer;


public class NodeIndex {

	public static class Node {

		final int statei;
		final int[] conf;
		final BigExp score;

		public Node(int statei, int[] conf, BigExp score) {
			this.statei = statei;
			this.conf = conf;
			this.score = score;
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Node && equals((Node)other);
		}

		public boolean equals(Node other) {
			return this.statei == other.statei
				&& Arrays.equals(this.conf, other.conf)
				&& this.score.equals(other.score);
		}

		@Override
		public String toString() {
			return String.format("Node[%d,%s,%s]", statei, Arrays.toString(conf), score);
		}
	}

	public final FixedDB db;
	public final String name;
	public final MultiStateConfSpace.State state;
	public final Consumer<Node> evictionListener;

	private final FixedIndex<BigExp,int[]> index;

	public NodeIndex(FixedDB db, String name, MultiStateConfSpace.State state, Consumer<Node> evictionListener) {

		this.db = db;
		this.name = name;
		this.state = state;
		this.evictionListener = evictionListener;

		index = new FixedIndex<>(
			db,
			name,
			Serializers.BigExp,
			Serializers.conf(state.confSpace),
			evictionListener == null ? null : (score, conf) ->
				evictionListener.accept(new Node(state.index, conf, score))
		);
	}

	public long size() {
		return index.size();
	}

	public BigExp lowestScore() {
		return index.lowestKey();
	}

	public BigExp highestScore() {
		return index.highestKey();
	}

	public void add(Node node) {
		assert (node.statei == state.index);
		index.add(node.score, node.conf);
	}

	public Node remove(BigExp score) {
		int[] conf = index.remove(score);
		return new Node(state.index, conf, score);
	}

	public void freeUpSpace() {
		index.freeUpSpace();
	}
}
