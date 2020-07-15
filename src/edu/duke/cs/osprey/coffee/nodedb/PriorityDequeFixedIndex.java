package edu.duke.cs.osprey.coffee.nodedb;

import com.google.common.collect.MinMaxPriorityQueue;

import java.nio.ByteBuffer;
import java.util.*;


/**
 * Trying to be a faster implementation of FixedIndex
 */
@SuppressWarnings("UnstableApiUsage") // the MinMaxPriorityQueue seems to be perpetually in "beta"... silly Google
public class PriorityDequeFixedIndex<S extends Comparable<S>, T extends FixedIndex.Indexable<S>> implements FixedIndex<S,T> {

	private static class Block<S extends Comparable<S>> {

		final long id;
		S max;

		Block(long id) {
			this.id = id;
			max = null;
		}

		<T extends Indexable<S>> void fill(BlockStore store, Serializer<T> serializer, MinMaxPriorityQueue<T> unpackedItems, int size) {

			max = null;

			// pack the lowest-scoring items
			ByteBuffer buf = store.get(id);
			buf.position(0);
			for (int i=0; i<size; i++) {
				var item = unpackedItems.pollFirst();
				if (i == size - 1) {
					max = item.score();
				}
				serializer.serialize(buf, item);
			}

			assert (max != null);
		}

		<T extends Indexable<S>> void moveAll(BlockStore store, Serializer<T> serializer, Collection<T> out, int size) {

			ByteBuffer buf = store.get(id);
			buf.position(0);
			for (int i=0; i<size; i++) {
				// TODO: this add is a huge bottleneck! 50% of the workload
				//  and half of that is BigExp comparisons
				checkAdd(out.add(serializer.deserialize(buf)));
			}
		}
	}

	public final BlockStore store;
	public final Serializer<T> serializer;

	public final int blockCapacity;
	public final Deque<T> dropped;

	// keep the blocks sorted by max val
	private final TreeSet<Block<S>> packedBlocks = new TreeSet<>(Comparator.comparing(block -> block.max));

	// keep the unpacked items sorted by score
	private final MinMaxPriorityQueue<T> unpackedItems;
	private final int unpackedCapacity;

	private long size = 0;

	public PriorityDequeFixedIndex(BlockStore store, Serializer<T> serializer) {

		this.store = store;
		this.serializer = serializer;

		blockCapacity = store.blockSize/serializer.bytes();
		dropped = new ArrayDeque<>(blockCapacity);
		unpackedCapacity = blockCapacity*2;
		unpackedItems = MinMaxPriorityQueue
			.orderedBy(Comparator.comparing((T item) -> item.score()))
			.maximumSize(unpackedCapacity + blockCapacity) // one extra block for temporary storage
			.create();
	}

	@Override
	public int blockCapacity() {
		return blockCapacity;
	}

	@Override
	public long size() {
		return size;
	}

	@Override
	public long freeSpace() {
		return store.numFreeBlocks()*blockCapacity + blockCapacity - unpackedItems.size();
	}

	private boolean packBlock() {

		// allocate a new block, if possible
		long blockid = store.allocateBlock();
		if (blockid == -1) {
			freeUpSpace();
			blockid = store.allocateBlock();
			if (blockid == -1) {
				return false;
			}
		}

		// fill it from the unpacked items and put it on the packed tree
		Block<S> block = new Block<>(blockid);
		block.fill(store, serializer, unpackedItems, blockCapacity);
		checkAdd(packedBlocks.add(block));
		return true;
	}

	@Override
	public boolean add(T item) {

		// if we hit capacity, pack a block
		if (unpackedItems.size() >= unpackedCapacity) {
			if (!packBlock()) {
				return false;
			}
		}

		// then add the item to the unpacked tree
		checkAdd(unpackedItems.add(item));
		size += 1;
		return true;
	}

	private boolean unpackHighestBlock() {

		// get the highest block
		Block<S> block = packedBlocks.pollLast();
		if (block == null) {
			return false;
		}

		block.moveAll(store, serializer, unpackedItems, blockCapacity);
		store.freeBlock(block.id);
		return true;
	}

	@Override
	public S highestScore() {

		S packedScore = null;
		if (!packedBlocks.isEmpty()) {
			packedScore = packedBlocks.last().max;
		}

		S unpackedScore = null;
		if (!unpackedItems.isEmpty()) {
			unpackedScore = unpackedItems.peekLast().score();
		}

		if (packedScore != null && unpackedScore != null) {

			if (packedScore.compareTo(unpackedScore) > 0) {
				return packedScore;
			} else {
				return unpackedScore;
			}

		} else if (packedScore != null) {
			return packedScore;
		} else {
			return unpackedScore;
		}
	}

	private T removeLastUnpacked() {
		var item = unpackedItems.pollLast();
		if (item == null) {
			throw new NoSuchElementException("out of unpacked items");
		}
		size -= 1;
		return item;
	}

	@Override
	public T removeHighest() {

		// if we're empty, there's nothing to remove
		if (unpackedItems.isEmpty() && packedBlocks.isEmpty()) {
			return null;
		}

		// if the highest item is already unpacked, return it
		if (packedBlocks.isEmpty() || (!unpackedItems.isEmpty() && unpackedItems.peekLast().score().compareTo(packedBlocks.last().max) > 0)) {
			return removeLastUnpacked();
		}

		// the highest item is packed, so unpack it first
		unpackHighestBlock();
		var item = removeLastUnpacked();

		if (unpackedItems.size() > unpackedCapacity) {
			// but don't overflow the unpacked tree
			boolean wasPacked = packBlock();
			assert (wasPacked);
		}

		return item;
	}

	@Override
	public void freeUpSpace() {

		// easy peasy, just remove the block with the lowest scores
		Block<S> block = packedBlocks.pollFirst();
		if (block == null) {
			return;
		}
		store.freeBlock(block.id);

		size -= blockCapacity;

		// put the items in the dropped queue
		if (dropped.isEmpty()) {
			block.moveAll(store, serializer, dropped, blockCapacity);
		} else {
			throw new IllegalStateException("dropped queue was not cleared before freeing up more space");
		}
	}

	@Override
	public Deque<T> dropped() {
		return dropped;
	}

	@Override
	public void clear() {

		for (var block : packedBlocks) {
			store.freeBlock(block.id);
		}
		packedBlocks.clear();
		size = 0;
	}

	private static void checkAdd(boolean wasAdded) {
		if (!wasAdded) {
			throw new IllegalStateException("item was not added");
		}
	}
}
