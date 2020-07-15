package edu.duke.cs.osprey.coffee.nodedb;

import com.hazelcast.cluster.Address;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;

import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;


/**
 * Tracks statistics about other members in the cluster.
 *
 * Totally thread-safe and suitable for a very highly multi-threaded environment.
 */
public class Neighbors {

	private class Neighbor {

		final Address addr;

		final long[] freeSpaces;
		final BigExp[] maxScores;
		long usedBytes;
		long totalBytes;

		Neighbor(Address addr) {
			this.addr = addr;
			freeSpaces = new long[confSpace.states.size()];
			maxScores = new BigExp[confSpace.states.size()];
		}

		<T> NeighborInfo<T> info(T item) {
			return new NeighborInfo<>(addr, item);
		}

		void receiveBroadcast(NodeIndices.BroadcastInfo nodeInfo) {

			// update the free spaces
			assert (nodeInfo.freeSpaces.length == confSpace.states.size());
			System.arraycopy(nodeInfo.freeSpaces, 0, this.freeSpaces, 0, confSpace.states.size());

			// update the max scores
			assert (nodeInfo.maxScores.length == confSpace.states.size());
			System.arraycopy(nodeInfo.maxScores, 0, this.maxScores, 0, confSpace.states.size());

			this.usedBytes = nodeInfo.usedBytes;
			this.totalBytes = nodeInfo.totalBytes;
		}

		void removeHighestNodes(int statei, int count, List<NodeIndex.Node> nodes) {
			NodeDB.checkSocketIOThread();
			var op = new GetHighestNodesOperation(statei, count);
			nodes.addAll(member.requestFrom(op, addr, 10, TimeUnit.SECONDS));
		}

		void addNodes(int statei, List<NodeIndex.Node> nodes) {
			NodeDB.checkSocketIOThread();
			var op = new AddNodesOperation(statei, nodes);
			member.sendTo(op, addr);
		}
	}

	public static class NeighborInfo<T> {

		public final Address addr;
		public final T item;

		public NeighborInfo(Address addr, T item) {
			this.addr = addr;
			this.item = item;
		}
	}

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;

	private final Map<Address,Neighbor> neighbors = new HashMap<>();

	public Neighbors(MultiStateConfSpace confSpace, ClusterMember member) {
		this.confSpace = confSpace;
		this.member = member;
	}

	private Neighbor getOrMake(Address addr) {
		return neighbors.computeIfAbsent(addr, key -> new Neighbor(addr));
	}

	public synchronized void receiveBroadcast(Address src, NodeIndices.BroadcastInfo nodeInfo) {
		getOrMake(src).receiveBroadcast(nodeInfo);
	}

	/**
	 * Returns the neighbor having the most free space,
	 * or null if no neighbors have any free space.
	 */
	public synchronized NeighborInfo<Long> findMostFreeSpace(int statei) {
		var neighbor = neighbors.values().stream()
			.filter(n -> n.freeSpaces[statei] > 0)
			.max(Comparator.comparing(n -> n.freeSpaces[statei]))
			.orElse(null);
		if (neighbor != null) {
			return neighbor.info(neighbor.freeSpaces[statei]);
		}
		return null;
	}

	/**
	 * Returns the neighbor having the highest nodes,
	 * or null if no neighbors have nodes.
	 */
	public synchronized NeighborInfo<BigExp> findHighestNodes(int statei) {
		var neighbor = neighbors.values().stream()
			.filter(n -> n.maxScores[statei] != null)
			.max(Comparator.comparing(n -> n.maxScores[statei]))
			.orElse(null);
		if (neighbor != null) {
			return neighbor.info(neighbor.maxScores[statei]);
		}
		return null;
	}

	public void addNodes(Address addr, int statei, List<NodeIndex.Node> nodes) {
		getOrMake(addr).addNodes(statei, nodes);
	}

	public void removeHighestNodes(Address addr, int statei, int count, List<NodeIndex.Node> nodes) {
		getOrMake(addr).removeHighestNodes(statei, count, nodes);
	}

	public float usage(long localUsedBytes, long localTotalBytes) {

		// start with the local usage
		long usedBytes = localUsedBytes;
		long totalBytes = localTotalBytes;

		// add usage from neighbors
		for (var n : neighbors.values()) {
			usedBytes += n.usedBytes;
			totalBytes += n.totalBytes;
		}

		// convert to a ratio
		return (float)usedBytes/(float)totalBytes;
	}
}
