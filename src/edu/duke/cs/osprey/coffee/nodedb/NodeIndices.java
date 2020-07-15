package edu.duke.cs.osprey.coffee.nodedb;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Streams;

import java.util.Arrays;
import java.util.List;
import java.util.function.Consumer;
import java.util.stream.Stream;


public class NodeIndices implements AutoCloseable {

	public static class BroadcastInfo {

		public final long[] freeSpaces;
		public final BigExp[] maxScores;
		public long usedBytes;
		public long totalBytes;

		public BroadcastInfo(int n) {
			freeSpaces = new long[n];
			maxScores = new BigExp[n];
		}

		public int size() {
			return freeSpaces.length;
		}
	}

	public final MultiStateConfSpace confSpace;
	public final long memBytes;

	public Consumer<Stream<NodeIndex.Node>> dropHandler = null;

	private final BlockStore store;
	private final NodeIndex[] indices;

	public NodeIndices(MultiStateConfSpace confSpace, long memBytes) {

		this.confSpace = confSpace;
		this.memBytes = memBytes;

		// allocate the block store
		store = new BlockStore(null, memBytes);

		// make sure there's at least 2 blocks for each index
		long minBytes = store.blockSize*confSpace.states.size()*2;
		if (memBytes < minBytes) {
			throw new IllegalArgumentException(String.format("NodeDB should have at least %d bytes for %d states",
				minBytes, confSpace.states.size()
			));
		}

		// init the state indices and metadata
		indices = confSpace.states.stream()
			.map(state -> new NodeIndex(store, state))
			.toArray(NodeIndex[]::new);
	}

	@Override
	public void close() {
		store.close();
	}

	public long size(int statei) {
		return indices[statei].size();
	}

	public void clear(int statei) {
		indices[statei].clear();
	}

	public BroadcastInfo getBroadcastInfo() {
		var out = new BroadcastInfo(indices.length);
		for (int i=0; i<indices.length; i++) {
			var index = indices[i];
			out.freeSpaces[i] = index.freeSpace();
			out.maxScores[i] = index.highestScore();
		}
		out.usedBytes = store.numUsedBytes();
		out.totalBytes = store.bytes;
		return out;
	}

	public void handleDropped() {

		// handle dropped nodes
		if (Arrays.stream(indices).anyMatch(index -> !index.dropped().isEmpty())) {

			// call the drop handler if possible
			if (dropHandler != null) {
				dropHandler.accept(Arrays.stream(indices)
					.flatMap(index -> index.dropped().stream())
				);
			}

			// forget the nodes
			for (var index : indices) {
				index.dropped().clear();
			}
		}
	}

	private void handleDropped(int statei) {

		var index = indices[statei];

		// any dropped nodes to process?
		if (index.dropped().isEmpty()) {
			return;
		}

		// call the drop handler if possible
		if (dropHandler != null) {
			dropHandler.accept(Streams.of(index.dropped()));
		}

		// then forget the nodes
		indices[statei].dropped().clear();
	}

	public BigExp highestScore(int statei) {
		return indices[statei].highestScore();
	}

	public boolean tryAdd(int statei, List<NodeIndex.Node> nodes) {

		// out of space?
		var index = indices[statei];
		if (nodes.size() > index.freeSpace()) {
			// yup
			return false;
		}

		// add the nodes
		for (var node : nodes) {
			boolean wasAdded = index.add(node);
			assert (wasAdded);
		}
		return true;
	}

	public void add(int statei, List<NodeIndex.Node> nodes) {

		var index = indices[statei];
		for (NodeIndex.Node node : nodes) {

			// add the node, if there's space
			boolean wasAdded = index.add(node);
			if (!wasAdded) {

				// free up space in all the other indices
				for (var state : confSpace.states) {
					if (state.index != node.statei) {
						indices[state.index].freeUpSpace();
					}
				}

				// try again
				wasAdded = indices[node.statei].add(node);
				if (!wasAdded) {
					throw new Error("Couldn't find/make space for a new node in the local store. This is a bug.");
				}
			}

			handleDropped();
		}
	}

	public void removeHighest(int statei, int count, List<NodeIndex.Node> nodes) {

		var index = indices[statei];

		for (int i=0; i<count; i++) {
			var node = index.removeHighest();

			// remove highest can drop nodes too, so handle them now
			handleDropped(statei);

			if (node == null) {
				break;
			}
			nodes.add(node);
		}
	}

	public long freeSpace(int statei) {
		return indices[statei].freeSpace();
	}

	public long nodesPerBlock(int statei) {
		return indices[statei].nodesPerBlock();
	}

	public long numUsedBytes() {
		return store.numUsedBytes();
	}
}
