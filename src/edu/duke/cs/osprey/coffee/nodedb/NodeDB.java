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

	public static final String ServiceName = "nodedb";
	private static final String ThreadName = "NodeDB";

	public static void checkSocketIOThread() {
		if (Thread.currentThread().getName().equals(NodeDB.ThreadName)) {
			throw new Error("don't do socket IO on the NodeDB thread");
		}
	}

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
	private final NodeIndices indices;
	private final RateLimitedThread broadcaster;
	private final Neighbors neighbors;

	private NodeDB(MultiStateConfSpace confSpace, ClusterMember member, File file, long fileBytes, long memBytes, long broadcastNs, File scoringLog) {

		this.confSpace = confSpace;
		this.member = member;
		this.file = file;
		this.fileBytes = fileBytes;
		this.memBytes = memBytes;
		this.broadcastNs = broadcastNs;
		this.scoringLog = scoringLog;

		// TODO: implement memory-buffered disk-backed options?
		// TEMP
		if (file != null || fileBytes > 0) {
			throw new Error("implement me");
		}

		perf = new NodePerformance(confSpace);
		perf.setLog(scoringLog);

		// the node indices aren't thread-safe, and can only be accessed by their creating thread
		// so make a thread to handle all the accesses
		thread = new BottleneckThread(ThreadName);
		indices = thread.get(() -> new NodeIndices(confSpace, memBytes));

		// make another thread to periodically keep the cluster members up-to-date
		broadcaster = new RateLimitedThread("NodeDB-bcast", broadcastNs, TimeUnit.NANOSECONDS, () -> broadcast());

		neighbors = new Neighbors(confSpace, member);

		// register NodeDB with hazelcast
		member.registerService(ServiceName, this);
		member.registerSerializer(NodeIndex.Node.class, Serializers.hazelcastNode(confSpace));

		// wait for everyone to catch up before sending the first broadcast
		member.barrier(1, TimeUnit.MINUTES);
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

		checkSocketIOThread();

		// get info from the indices
		var info = thread.get(() -> indices.getBroadcastInfo());

		// broadcast
		member.sendToOthers(() -> new BroadcastOperation(info, perf));
	}

	void receiveBroadcast(Address src, NodeIndices.BroadcastInfo nodeInfo) {
		neighbors.receiveBroadcast(src, nodeInfo);
	}

	/**
	 * Remove all the nodes from the given state
	 */
	public void clear(int statei) {

		// propagate to neighbors
		checkSocketIOThread();
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
		var neighbor = neighbors.findMostFreeSpace(statei);
		if (neighbor != null) {
			neighbors.addNodes(neighbor.addr, statei, nodes);
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
		var neighbor = thread.get(() -> {

			// compare the local scores with the highest neighbor to figure out where the best nodes are
			var highestNeighbor = neighbors.findHighestNodes(statei);
			BigExp localMaxScore = indices.highestScore(statei);

			// if the local nodes the best, get those
			if (localMaxScore != null && (highestNeighbor == null || localMaxScore.compareTo(highestNeighbor.item) > 0)) {
				indices.removeHighest(statei, count, nodes);
				broadcaster.request();
				return null;
			}

			return highestNeighbor;
		});

		// if the neighbor's nodes are the best, get those
		if (neighbor != null) {
			neighbors.removeHighestNodes(neighbor.addr, statei, count, nodes);
		}

		// otherwise, there are no nodes anywhere
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
		return neighbors.usage(
			thread.get(() -> indices.numUsedBytes()),
			memBytes
		);
	}
}
