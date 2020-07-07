package edu.duke.cs.osprey.coffee.nodedb;

import java.nio.ByteBuffer;
import java.util.ArrayDeque;
import java.util.Comparator;
import java.util.Deque;
import java.util.TreeSet;


/**
 * A very simple implementation that turns out to have slow removeHighest() speeds. =(
 * Let's just call it a "reference" implentation and move on.
 */
public class ReferenceFixedIndex<S extends Comparable<S>, T extends FixedIndex.Indexable<S>> implements FixedIndex<S,T> {

	private static class Block<S extends Comparable<S>> {

		final long id;
		int size;
		S min;
		S max;

		public Block(long id) {
			this.id = id;
			size = 0;
			min = null;
			max = null;
		}

		public <T extends FixedIndex.Indexable<S>> void add(BlockStore store, FixedIndex.Serializer<T> serializer, T item) {

			// append the item to the block
			ByteBuffer buf = store.get(id);
			int bytes = serializer.bytes();
			buf.position(bytes*size);
			serializer.serialize(buf, item);

			// update block metadata
			size += 1;
			S score = item.score();
			if (min == null || score.compareTo(min) < 0) {
				min = score;
			}
			if (max == null || score.compareTo(max) > 0) {
				max = score;
			}
		}

		public <T extends FixedIndex.Indexable<S>> T removeHighest(BlockStore store, FixedIndex.Serializer<T> serializer) {

			ByteBuffer buf = store.get(id);

			T bestItem = null;
			S secondBestScore = null;
			int bestIndex = 0;

			// find the best and second-best items in the block
			buf.position(0);
			for (int i=0; i<size; i++) {
				T item = serializer.deserialize(buf);
				if (bestItem == null) {
					bestItem = item;
					bestIndex = i;
				} else if (item.score().compareTo(bestItem.score()) > 0) {
					secondBestScore = bestItem.score();
					bestItem = item;
					bestIndex = i;
				} else if (secondBestScore == null || item.score().compareTo(secondBestScore) > 0) {
					secondBestScore = item.score();
				}
			}

			// remove it, shift the remaining items down
			int itemSize = serializer.bytes();
			for (int i=bestIndex; i<size - 1; i++) {
				buf.position((i + 1)*itemSize);
				T item = serializer.deserialize(buf);
				buf.position(i*itemSize);
				serializer.serialize(buf, item);
			}
			size -= 1;
			max = secondBestScore;

			return bestItem;
		}

		public <T extends FixedIndex.Indexable<S>> void getAll(BlockStore store, FixedIndex.Serializer<T> serializer, Deque<T> dropped) {

			ByteBuffer buf = store.get(id);

			buf.position(0);
			for (int i=0; i<size; i++) {
				dropped.add(serializer.deserialize(buf));
			}
		}
	}

	public final BlockStore store;
	public final FixedIndex.Serializer<T> serializer;

	public final int blockCapacity;
	public final Deque<T> dropped;

	// keep the blocks sorted by max val
	private final TreeSet<Block<S>> packedBlocks = new TreeSet<>(Comparator.comparing(block -> block.max));

	private Block<S> unpackedBlock = null;
	private long size = 0;

	public ReferenceFixedIndex(BlockStore store, FixedIndex.Serializer<T> serializer) {

		this.store = store;
		this.serializer = serializer;

		blockCapacity = store.blockSize/serializer.bytes();
		dropped = new ArrayDeque<>(blockCapacity);
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
		long freeSpace = store.numFreeBlocks()*blockCapacity;
		if (unpackedBlock != null) {
			freeSpace += blockCapacity - unpackedBlock.size;
		}
		return freeSpace;
	}

	@Override
	public boolean add(T item) {

		if (unpackedBlock == null) {

			// allocate a new block
			long blockid = store.allocateBlock();
			if (blockid == -1) {
				freeUpSpace();
				blockid = store.allocateBlock();
				if (blockid == -1) {
					return false;
				}
			}
			unpackedBlock = new Block<>(blockid);
		}

		assert (unpackedBlock.size < blockCapacity);

		// add the item to the block
		unpackedBlock.add(store, serializer, item);

		assert (unpackedBlock.size <= blockCapacity);

		// if we packed the block, put it on the tree
		if (unpackedBlock.size == blockCapacity) {
			packedBlocks.add(unpackedBlock);
			unpackedBlock = null;
		}

		size += 1;

		return true;
	}

	@Override
	public S highestScore() {

		S packedScore = null;
		if (!packedBlocks.isEmpty()) {
			packedScore = packedBlocks.last().max;
		}

		S unpackedScore = null;
		if (unpackedBlock != null) {
			unpackedScore = unpackedBlock.max;
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

	@Override
	public T removeHighest() {

		if (unpackedBlock == null) {

			// find the block with the highest scores
			Block<S> block = packedBlocks.pollLast();
			if (block == null) {
				return null;
			}

			// remove the highest-scoring item from the block
			T item = block.removeHighest(store, serializer);
			size -= 1;
			assert (block.size > 0);
			unpackedBlock = block;

			return item;

		} else {

			// find the block with the highest scores
			if (!packedBlocks.isEmpty() && packedBlocks.last().max.compareTo(unpackedBlock.max) > 0) {

				// it's a packed block
				Block<S> block = packedBlocks.pollLast();

				// replace its highest-scoring item with something from the unpacked block
				T item = block.removeHighest(store, serializer);
				size -= 1;
				block.add(store, serializer, unpackedBlock.removeHighest(store, serializer));
				assert (block.size == blockCapacity);
				if (unpackedBlock.size == 0) {
					store.freeBlock(unpackedBlock.id);
					unpackedBlock = null;
				}
				packedBlocks.add(block);

				return item;

			} else {

				// it's the unpacked block, remove its highest item
				T item = unpackedBlock.removeHighest(store, serializer);
				size -= 1;
				if (unpackedBlock.size == 0) {
					store.freeBlock(unpackedBlock.id);
					unpackedBlock = null;
				}

				return item;
			}
		}
	}

	@Override
	public void freeUpSpace() {

		// easy peasy, just remove the block with the lowest scores
		Block<S> block = packedBlocks.pollFirst();
		store.freeBlock(block.id);

		size -= block.size;

		// put the items in the dropped queue
		if (dropped.isEmpty()) {
			block.getAll(store, serializer, dropped);
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
}
