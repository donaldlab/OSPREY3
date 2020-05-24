package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.tools.BigExp;
import org.mapdb.*;

import java.io.File;
import java.util.Arrays;


public class NodeDB {

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;
	public final File file;

	public NodeDB(MultiStateConfSpace confSpace, ClusterMember member, File file) {

		this.confSpace = confSpace;
		this.member = member;
		this.file = file;

		// TODO
	}

	public void init(int statei, int[] conf, BigExp score) {

		// TODO: update the index and the node store

		// storage
		// use MapDB HTreeMap: https://jankotek.gitbooks.io/mapdb/htreemap/
		// allocate X size to an in-memory index,
			// use a BTree map with a fixed-size volume: ByteBufferMemoryVolSingle
				// test behavior when running out of space!!
			// just need to keep track of confs with big-ish zSumUpper values, but not a whole sorted list
			// and do the same for small zSumUpper values, for eviction
			// if in-memory index is exhausted, rebuild it by scanning the hashmap

		/* TEMP
		DB dbFile = DBMaker
			.fileDB(new File(""))
			.make();
		var mapFile = dbFile.hashMap("foo")
			.keySerializer(Serializer.STRING)
			.valueSerializer(Serializer.STRING)
			.expireStoreSize(5)
			.expireAfterUpdate()
			.open();
		*/

		/* TEMP
		// one index per state
		DB dbIndex = DBMaker
			.volumeDB(null, false) // TODO: copy ByteBufferVolSingle and ByteBufferMemoryVolSingle
			.make();
		var index = dbIndex.treeMap("index")
			.open();

		StoreDirect store = (StoreDirect)index.getStore();
		long freeBytes = store.getFreeSize();
		long bytes = store.getTotalSize();
		*/

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
}
