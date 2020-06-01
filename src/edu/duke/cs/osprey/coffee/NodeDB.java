package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import org.mapdb.CC;
import org.mapdb.DBException;

import java.io.File;


public class NodeDB {

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;
	public final File file;
	public final long fileBytes;
	public final long memBytes;

	private final FixedDB db;
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

		// make sure there's at least 1 page for the db, and 2 pages for each index
		long minBytes = CC.PAGE_SIZE*(1 + 2*confSpace.states.size());
		if (memBytes < minBytes) {
			throw new IllegalArgumentException(String.format("NodeDB should have at least %d bytes for %d states",
				minBytes, confSpace.states.size()
			));
		}

		db = new FixedDB(null, memBytes);
		indices = confSpace.states.stream()
			.map(state -> new NodeIndex(db, "nodeindex-" + state.name, state, null))
			.toArray(NodeIndex[]::new);
	}

	public long size(int statei) {
		return indices[statei].size();
	}

	/**
	 * Add the node to the local store
	 */
	public void addLocal(NodeIndex.Node node) {
		try {

			// add the node, if there's space
			indices[node.statei].add(node);

		} catch (DBException.VolumeMaxSizeExceeded ex) {

			// free up space in all the other indices
			for (var state : confSpace.states) {
				if (state.index != node.statei) {
					indices[state.index].freeUpSpace();
				}
			}

			try {
				// try again
				indices[node.statei].add(node);

			} catch (DBException.VolumeMaxSizeExceeded ex2) {
				throw new Error("Couldn't find/make space for a new node in the local store. This is a bug.");
			}
		}
	}

	public NodeIndex.Node pollHighestLocal(int statei) {
		var index = indices[statei];
		return index.remove(index.highestScore());
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
