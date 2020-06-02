package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.db.BlockStore;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;

import java.io.File;


public class NodeDB {

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;
	public final File file;
	public final long fileBytes;
	public final long memBytes;

	private final BlockStore store;
	private final NodeIndex[] indices;

	public NodeDB(MultiStateConfSpace confSpace, ClusterMember member, File file, long fileBytes, long memBytes) {

		this.confSpace = confSpace;
		this.member = member;
		this.file = file;
		this.fileBytes = fileBytes;
		this.memBytes = memBytes;

		// TODO: implement memory-buffered disk-backed options
		// TEMP
		if (file != null || fileBytes > 0) {
			throw new Error("implement me");
		}

		// allocate the block store
		store = new BlockStore(null, memBytes);

		// make sure there's at least 2 blocks for each index
		long minBytes = store.blockSize*confSpace.states.size()*2;
		if (memBytes < minBytes) {
			throw new IllegalArgumentException(String.format("NodeDB should have at least %d bytes for %d states",
				minBytes, confSpace.states.size()
			));
		}

		indices = confSpace.states.stream()
			.map(state -> new NodeIndex(store, state))
			.toArray(NodeIndex[]::new);
	}

	public long size(int statei) {
		return indices[statei].size();
	}

	/**
	 * Add the node to the local store
	 */
	public void addLocal(NodeIndex.Node node) {

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
	}

	public NodeIndex.Node removeHighestLocal(int statei) {
		var index = indices[statei];
		return index.removeHighest();
	}

	// querying node:
		// query node from random member that has good nodes
			// in proportion to the goodness of the nodes
			// requires everyone to track everyone else's goodness
		// send copy of node to requesting member
		// mark node as processing, exclude it from future queries

	// saving nodes:
		// if there's space, just save locally
		// if not, pick another member to save them
			// pick a member with free space, at random
			// if no one has free space, pick any random member
		// that member must save the nodes, at least the good nodes
		// that member must evict nodes until there's enough space
		// requires every member to track every other member's free space

	// node eviction:
		// evict nodes with low zSumUpper values
			// relative to the root zSumUpper for the state
		// rely on index to find them
}
