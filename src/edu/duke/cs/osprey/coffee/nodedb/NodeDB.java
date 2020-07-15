package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.cluster.Address;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.coffee.Serializers;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.parallelism.BottleneckThread;
import edu.duke.cs.osprey.parallelism.RateLimitedThread;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.File;
import java.util.*;
import java.util.concurrent.TimeUnit;
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
		private File scoringLog;

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

		public Builder setScoringLog(File val) {
			scoringLog = val;
			return this;
		}

		public NodeDB build() {
			return new NodeDB(
				confSpace,
				member,
				file, fileBytes,
				memBytes,
				broadcastNs,
				scoringLog
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

		public void removeHighest(int statei, int count, List<NodeIndex.Node> nodes) {
			nodes.addAll(member.requestFrom(new GetHighestNodesOperation(statei, count), addr, 10, TimeUnit.SECONDS));
		}

		public void add(int statei, List<NodeIndex.Node> nodes) {
			member.sendTo(new AddNodesOperation(statei, nodes), addr);
		}
	}

	public static final String ServiceName = "nodedb";

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;
	public final File file;
	public final long fileBytes;
	public final long memBytes;
	public final long broadcastNs;
	public final File scoringLog;

	public final NodePerformance perf;

	/**
	 * Since the NodeIndex instances are not thread-safe,
	 * and the BlockStore memory must be accessed by a single thread,
	 * we have to serialize all DB accesses through one thread.
	 *
	 * This could eventually end up being too slow,
	 * but profiling shows the performance isn't too bad yet,
	 * even on 48 threads.
	 */
	private final BottleneckThread thread;
	private final RateLimitedThread broadcaster;

	private NodeIndices indices;
	private Map<Address,Neighbor> neighbors;

	private NodeDB(MultiStateConfSpace confSpace, ClusterMember member, File file, long fileBytes, long memBytes, long broadcastNs, File scoringLog) {

		this.confSpace = confSpace;
		this.member = member;
		this.file = file;
		this.fileBytes = fileBytes;
		this.memBytes = memBytes;
		this.broadcastNs = broadcastNs;
		this.scoringLog = scoringLog;

		perf = new NodePerformance(confSpace);
		perf.setLog(scoringLog);

		// TODO: implement memory-buffered disk-backed options
		// TEMP
		if (file != null || fileBytes > 0) {
			throw new Error("implement me");
		}

		// the node indices aren't thread-safe, and can only be accessed by their creating thread
		// so make a thread to handle all the accesses
		thread = new BottleneckThread("NodeDB");

		// make another thread to periodically keep the cluster members up-to-date
		broadcaster = new RateLimitedThread("NodeDB-bcast", broadcastNs, TimeUnit.NANOSECONDS, () -> broadcast());

		thread.exec(() -> {

			indices = new NodeIndices(confSpace, memBytes);

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
		});

		// wait for everyone to catch up before sending the first broadcast
		member.barrier(1, TimeUnit.MINUTES);

		// but tell them that we care
		broadcast();
	}

	/**
	 * Set a function to call when dropped nodes need to be processed.
	 * Called from the NodeDB thread, not the caller thread!
	 **/
	public void setDropHandler(Consumer<Stream<NodeIndex.Node>> dropHandler) {
		thread.exec(() -> indices.dropHandler = dropHandler);
	}

	@Override
	public void close() {
		broadcaster.close();
		thread.exec(() -> indices.close());
		thread.close();
	}

	public long size(int statei) {
		return thread.get(() -> indices.size(statei));
	}

	public void broadcast() {

		// get info from the indices
		var info = thread.get(() -> indices.getBroadcastInfo());

		// broadcast
		// NOTE: don't send on the db thread, we don't want that thread blocking on network IO
		member.sendToOthers(() -> new BroadcastOperation(info, perf));
	}

	void receiveBroadcast(Address src, long[] freeSpaces, BigExp[] maxScores, long usedBytes, long totalBytes) {
		// TODO: do we need to synchronize neighbors?
		getNeighbor(src).broadcasted(freeSpaces, maxScores, usedBytes, totalBytes);
	}

	private Neighbor getNeighbor(Address addr) {
		return neighbors.computeIfAbsent(addr, key -> new Neighbor(addr));
	}

	/**
	 * Remove all the nodes from the given state
	 */
	public void clear(int statei) {

		// propagate to neighbors
		member.sendToOthers(() -> new ClearOperation(statei));

		clearLocal(statei);
	}

	public void clearLocal(int statei) {
		thread.exec(() -> {
			indices.clear(statei);
			broadcaster.request();
		});
	}

	/**
	 * Adds nodes to the cluster.
	 * Local storage is preferred if there's space.
	 * Next, remote storage is preferred if there's space.
	 * Otherwise, space will be evicted from local storage to make room.
	 */
	public void add(int statei, List<NodeIndex.Node> nodes) {

		// prefer local storage first
		boolean wasAdded = thread.get(() -> indices.tryAdd(statei, nodes));
		if (wasAdded) {
			broadcaster.request();
			return;
		}

		// prefer remote storage next
		Neighbor neighbor = neighbors.values().stream()
			.max(Comparator.comparing(n -> n.freeSpaces[statei]))
			.orElse(null);
		if (neighbor != null && neighbor.freeSpaces[statei] > 0) {
			neighbor.add(statei, nodes);
			return;
		}

		// finally, force local storage
		thread.exec(() -> indices.add(statei, nodes));
		broadcaster.request();
	}

	/**
	 * Conveience method to add a single node.
	 * The batched version is preferred, for speed.
	 */
	public void add(NodeIndex.Node node) {
		add(node.statei, Collections.singletonList(node));
	}

	/**
	 * Add nodes to the local store
	 */
	public void addLocal(int statei, List<NodeIndex.Node> nodes) {
		thread.exec(() -> {
			indices.add(statei, nodes);
			broadcaster.request();
		});
	}

	/**
	 * Conveience method to add a single node to the local store.
	 * The batched version is preferred, for speed.
	 */
	public void addLocal(NodeIndex.Node node) {
		addLocal(node.statei, Collections.singletonList(node));
	}

	/**
	 * Removes the highest node from the local index.
	 */
	public void removeHighestLocal(int statei, int count, List<NodeIndex.Node> nodes) {
		thread.exec(() -> {
			indices.removeHighest(statei, count, nodes);
			broadcaster.request();
		});
	}

	/**
	 * Conveience method to remove a single node.
	 * The batched version is preferred, for speed.
	 */
	public NodeIndex.Node removeHighestLocal(int statei) {
		var nodes = new ArrayList<NodeIndex.Node>(1);
		removeHighestLocal(statei, 1, nodes);
		return nodes.get(0);
	}

	/**
	 * Removes high-scoring nodes from the cluster.
	 *
	 * It's not necessarily the highest-scoring nodes,
	 * due to delays in score broadcasting,
	 * but they should be pretty high.
	 */
	public void removeHigh(int statei, int count, List<NodeIndex.Node> nodes) {
		thread.exec(() -> {

			// find the neighbor with the highest node for this state
			Neighbor highestNeighbor = neighbors.values().stream()
				.filter(neighbor -> neighbor.maxScores[statei] != null)
				.max(Comparator.comparing(neighbor -> neighbor.maxScores[statei]))
				.orElse(null);

			BigExp localMaxScore = indices.highestScore(statei);

			if (highestNeighbor != null) {

				// is a local node higher?
				if (localMaxScore != null && localMaxScore.compareTo(highestNeighbor.maxScores[statei]) > 0) {

					// yup, use that
					indices.removeHighest(statei, count, nodes);
					broadcaster.request();

				} else {

					// nope, use the cluster node
					highestNeighbor.removeHighest(statei, count, nodes);
				}

			} else if (localMaxScore != null) {

				// use the local node
				indices.removeHighest(statei, count, nodes);
				broadcaster.request();

			//} else {

				// no nodes anywhere
			}
		});
	}

	/**
	 * Conveience method to remove a single node.
	 * The batched version is preferred, for speed.
	 */
	public NodeIndex.Node removeHigh(int statei) {
		var nodes = new ArrayList<NodeIndex.Node>(1);
		removeHigh(statei, 1, nodes);
		return nodes.get(0);
	}

	public long freeSpaceLocal(int statei) {
		return thread.get(() -> indices.freeSpace(statei));
	}

	public long nodesPerBlock(int statei) {
		// constant lookup, don't need to synchronize
		return indices.nodesPerBlock(statei);
	}

	/**
	 * Returns the ratio of used space to total space.
	 */
	public float usage() {

		// start with the local usage
		long usedBytes = thread.get(() -> indices.numUsedBytes());
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
