package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.db.BlockStore;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.IntEncoding;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.stream.IntStream;


public class NodeIndex {

	public static class Node implements FixedIndex.Indexable<BigExp> {

		final int statei;
		final int[] conf;
		final BigExp score;

		public Node(int statei, int[] conf, BigExp score) {
			this.statei = statei;
			this.conf = conf;
			this.score = score;
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
				&& this.score.equals(other.score);
		}

		@Override
		public String toString() {
			return String.format("Node[%d,%s,%s]", statei, Arrays.toString(conf), score);
		}

		private static class Serializer implements FixedIndex.Serializer<Node> {

			final MultiStateConfSpace.State state;
			final IntEncoding confEncoding;

			Serializer(MultiStateConfSpace.State state) {
				this.state = state;
				confEncoding = IntEncoding.get(
					IntStream.range(0, state.confSpace.numPos())
						.map(posi -> state.confSpace.numConf(posi) - 1)
						.max()
						.orElse(-1)
				);
			}

			@Override
			public int bytes() {
				return 4+8 // BigExp
					+ state.confSpace.numPos()*confEncoding.numBytes; // conf
			}

			@Override
			public void serialize(ByteBuffer out, Node node) {
				out.putInt(node.score.exp);
				out.putDouble(node.score.fp);
				for (int i=0; i<state.confSpace.numPos(); i++) {
					switch (confEncoding) {
						case Byte -> out.put((byte)node.conf[i]);
						case Short -> out.putShort((short)node.conf[i]);
						case Int -> out.putInt(node.conf[i]);
					}
				}
			}

			@Override
			public Node deserialize(ByteBuffer in) {
				int exp = in.getInt();
				double fp = in.getDouble();
				BigExp score = new BigExp(fp, exp);
				int[] conf = new int[state.confSpace.numPos()];
				for (int i=0; i<conf.length; i++) {
					switch (confEncoding) {
						case Byte -> conf[i] = in.get();
						case Short -> conf[i] = in.getShort();
						case Int -> conf[i] = in.getInt();
					}
				}
				return new Node(state.index, conf, score);
			}
		}
	}

	public final BlockStore store;
	public final MultiStateConfSpace.State state;

	private final FixedIndex<BigExp,Node> index;

	public NodeIndex(BlockStore store, MultiStateConfSpace.State state) {

		this.store = store;
		this.state = state;

		index = new FixedIndex<>(store, new Node.Serializer(state));
	}

	protected int nodesPerBlock() {
		return index.blockCapacity;
	}

	public long size() {
		return index.size();
	}

	public boolean add(Node node) {
		assert (node.statei == state.index);
		return index.add(node);
	}

	public Node removeHighest() {
		return index.removeHighest();
	}

	public void freeUpSpace() {
		index.freeUpSpace();
	}
}
