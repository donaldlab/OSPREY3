package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.cluster.Address;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.File;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class NodeDB implements AutoCloseable {

	// oh, how I wish Java had defaultable arguments ...

	public static class Builder {

		public final MultiStateConfSpace confSpace;
		public final ClusterMember member;

		private File file = null;
		private long fileBytes = 0;
		private long memBytes = 0;
		private long broadcastNs = 1_000_000_000L; // 1 second

		public Builder(MultiStateConfSpace confSpace, ClusterMember member) {
			this.confSpace = confSpace;
			this.member = member;
		}

		public Builder setFile(File file, long bytes) {
			this.file = file;
			this.fileBytes = bytes;
			return this;
		}

		public Builder setMem(long bytes) {
			this.memBytes = bytes;
			return this;
		}

		public Builder setBroadcastNs(long ns) {
			this.broadcastNs = ns;
			return this;
		}

		public NodeDB build() {
			return new NodeDB(
				confSpace,
				member,
				file, fileBytes,
				memBytes,
				broadcastNs
			);
		}
	}

	private class Neighbor {

		final Address addr;
		final long[] freeSpaces;
		final BigExp[] maxScores;

		long usedBytes;
		long totalBytes;

		public Neighbor(Address addr) {
			this.addr = addr;
			freeSpaces = new long[confSpace.states.size()];
			maxScores = new BigExp[confSpace.states.size()];
		}

		public void broadcasted(long[] freeSpaces, BigExp[] maxScores, long usedBytes, long totalBytes) {

			// update the free spaces
			assert (freeSpaces.length == confSpace.states.size());
			System.arraycopy(freeSpaces, 0, this.freeSpaces, 0, confSpace.states.size());

			// update the max scores
			assert (maxScores.length == confSpace.states.size());
			System.arraycopy(maxScores, 0, this.maxScores, 0, confSpace.states.size());

			this.usedBytes = usedBytes;
			this.totalBytes = totalBytes;
		}

		public NodeIndex.Node removeHighest(int statei) {
			return member.requestFrom(new GetHighestNodeOperation(statei), addr, 10, TimeUnit.SECONDS);
		}

		public void add(NodeIndex.Node node) {
			member.sendTo(new AddNodeOperation(node), addr);
		}
	}

	public static final String ServiceName = "nodedb";

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;
	public final File file;
	public final long fileBytes;
	public final long memBytes;
	public final long broadcastNs;

	/**
	 * Function to call when dropped nodes need to be processed.
	 * Called from the NodeDB thread, not the caller thread!
	 **/
	public Consumer<Stream<NodeIndex.Node>> dropHandler = null;

	public final NodePerformance perf;

	private final ThreadPoolTaskExecutor tasks;
	private Map<Address,Neighbor> neighbors;
	private BlockStore store;
	private NodeIndex[] indices;

	private long lastBroadcastNs = 0;

	private void threadExec(Runnable task) {
		tasks.submit(
			() -> {
				task.run();
				return 42;
			},
			answer -> {}
		);
		tasks.waitForFinish();
	}

	private <T> T threadGet(TaskExecutor.Task<T> task) {
		var ref = new AtomicReference<T>(null);
		tasks.submit(task, ref::set);
		tasks.waitForFinish();
		return ref.get();
	}

	private NodeDB(MultiStateConfSpace confSpace, ClusterMember member, File file, long fileBytes, long memBytes, long broadcastNs) {

		this.confSpace = confSpace;
		this.member = member;
		this.file = file;
		this.fileBytes = fileBytes;
		this.memBytes = memBytes;
		this.broadcastNs = broadcastNs;

		perf = new NodePerformance(confSpace);

		// TODO: implement memory-buffered disk-backed options
		// TEMP
		if (file != null || fileBytes > 0) {
			throw new Error("implement me");
		}

		// the node indices aren't thread-safe, and can only be accessed by their creating thread
		// so make a thread to handle all the accesses
		tasks = new ThreadPoolTaskExecutor();
		tasks.start(1);

		threadExec(() -> {

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

			// register NodeDB with hazelcast
			member.registerService(ServiceName, this);
			member.registerSerializer(NodeIndex.Node.class, Serializers.hazelcastNode(confSpace));

			// wait for everyone to finish registering with hazelcast, so NodeDB operations will work
			member.barrier(1, TimeUnit.MINUTES);

			// spy on our neighbors
			neighbors = member.otherMemberAddresses().stream()
				.collect(Collectors.toMap(
					addr -> addr,
					addr -> new Neighbor(addr)
				));

			// but tell them that we care
			broadcast();
		});
	}

	@Override
	public void close() {
		threadExec(() ->
			store.close()
		);
		tasks.close();
	}

	public long size(int statei) {
		return indices[statei].size();
	}

	/**
	 * Add the node to the local store
	 */
	public void addLocal(NodeIndex.Node node) {
		threadExec(() -> addLocalOnThread(node));
	}

	private void addLocalOnThread(NodeIndex.Node node) {

		// add the node, if there's space
		boolean wasAdded = indices[node.statei].add(node);
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

		broadcastIfNeeded();
	}

	private void broadcastIfNeeded() {

		// is it time to check again?
		if (System.nanoTime() - lastBroadcastNs < broadcastNs) {

			// nope
			return;
		}

		broadcast();
	}

	public void broadcast() {

		// get the free spaces
		long[] freeSpaces = Arrays.stream(indices)
			.mapToLong(NodeIndex::freeSpace)
			.toArray();

		// get the max scores
		BigExp[] maxScores = Arrays.stream(indices)
			.map(NodeIndex::highestScore)
			.toArray(BigExp[]::new);

		// broadcast
		member.sendToOthers(() -> new BroadcastOperation(freeSpaces, maxScores, store.numUsedBytes(), memBytes, perf));
		lastBroadcastNs = System.nanoTime();
	}

	public void receiveBroadcast(Address src, long[] freeSpaces, BigExp[] maxScores, long usedBytes, long totalBytes) {
		threadExec(() ->
			getNeighbor(src).broadcasted(freeSpaces, maxScores, usedBytes, totalBytes)
		);
	}

	private Neighbor getNeighbor(Address addr) {
		return neighbors.computeIfAbsent(addr, key -> new Neighbor(addr));
	}

	public NodeIndex.Node removeHighestLocal(int statei) {
		return threadGet(() -> {

			var index = indices[statei];
			NodeIndex.Node node = index.removeHighest();

			broadcastIfNeeded();

			return node;
		});
	}

	/**
	 * Removes a high-scoring node from the cluster.
	 *
	 * It's not necessarily the highest-scoring node,
	 * due to delays in score broadcasting,
	 * but it should be pretty high.
	 */
	public NodeIndex.Node removeHigh(int statei) {
		return threadGet(() -> {

			// find the neighbor with the highest node for this state
			Neighbor highestNeighbor = neighbors.values().stream()
				.filter(neighbor -> neighbor.maxScores[statei] != null)
				.max(Comparator.comparing(neighbor -> neighbor.maxScores[statei]))
				.orElse(null);

			BigExp localMaxScore = indices[statei].highestScore();

			if (highestNeighbor != null) {

				// is a local node higher?
				if (localMaxScore != null && localMaxScore.compareTo(highestNeighbor.maxScores[statei]) > 0) {

					// yup, use that
					return indices[statei].removeHighest();

				} else {

					// nope, use the cluster node
					return highestNeighbor.removeHighest(statei);
				}

			} else if (localMaxScore != null) {

				// use the local node
				return indices[statei].removeHighest();

			} else {

				// no nodes anywhere
				return null;
			}
		});
	}

	/**
	 * Adds the node to the cluster.
	 * Local storage is preferred if there's space.
	 * Next, remote storage is preferred if there's space.
	 * Otherwise, space will be evicted from local storage to make room.
	 */
	public void add(NodeIndex.Node node) {
		threadExec(() -> {

			// prefer local storage first
			if (indices[node.statei].freeSpace() > 0) {
				addLocalOnThread(node);
				return;
			}

			// prefer remote storage next
			Neighbor neighbor = neighbors.values().stream()
				.max(Comparator.comparing(n -> n.freeSpaces[node.statei]))
				.orElse(null);
			if (neighbor != null && neighbor.freeSpaces[node.statei] > 0) {
				neighbor.add(node);
				return;
			}

			// finally, force local storage
			addLocalOnThread(node);
		});
	}

	public long freeSpaceLocal(int statei) {
		return threadGet(() ->
			indices[statei].freeSpace()
		);
	}

	public long nodesPerBlock(int statei) {
		// constant lookup, don't need to synchronize
		return indices[statei].nodesPerBlock();
	}

	/**
	 * Returns the ratio of used space to total space.
	 */
	public float usage() {

		// start with the local usage
		long usedBytes = store.numUsedBytes();
		long totalBytes = memBytes;

		// add usage from neighbors
		for (var n : neighbors.values()) {
			usedBytes += n.usedBytes;
			totalBytes += n.totalBytes;
		}

		// convert to a ratio
		return (float)usedBytes/(float)totalBytes;
	}
}
