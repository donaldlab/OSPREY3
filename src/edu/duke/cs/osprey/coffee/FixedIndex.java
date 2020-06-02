package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.db.BlockStore;

import java.nio.ByteBuffer;
import java.util.*;


/**
 * An index that allows queries for items with high scores,
 * but discards items with low scores when it runs out of space.
 *
 * Performance-wise, the index is designed for fast add (and freeUpSpace) speeds,
 * perhaps at the expense of remove speeds.
 *
 * This implementation is NOT thread-safe!
 */
public class FixedIndex<S extends Comparable<S>, T extends FixedIndex.Indexable<S>> {

	// TODO: make thread-safe?

	public interface Indexable<S> {
		S score();
	}

	public interface Serializer<T> {
		int bytes();
		void serialize(ByteBuffer out, T data);
		T deserialize(ByteBuffer in);
	}

	private static class Block<S extends Comparable<S>> {

		final long id;
		int size;
		S min;
		S max;

		public Block(long id) {
			this.id = id;
			size = 1;
			min = null;
			max = null;
		}

		public <T extends Indexable<S>> void add(BlockStore store, Serializer<T> serializer, T item) {

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

		public <T extends Indexable<S>> T removeHighest(BlockStore store, Serializer<T> serializer) {

			ByteBuffer buf = store.get(id);
			int itemSize = serializer.bytes();

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
	}

	public final BlockStore store;
	public final Serializer<T> serializer;

	public final int blockCapacity;

	// keep the blocks sorted by max val
	private final TreeSet<Block<S>> packedBlocks = new TreeSet<>(Comparator.comparing(block -> block.max));

	private Block<S> unpackedBlock = null;
	private long size = 0;

	public FixedIndex(BlockStore store, Serializer<T> serializer) {

		this.store = store;
		this.serializer = serializer;

		blockCapacity = store.blockSize/serializer.bytes();
	}

	public long size() {
		return size;
	}

	/**
	 * Adds the item to the index.
	 * Returns true if the item was added, false if there wasn't enough space.
	 */
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

	public void freeUpSpace() {

		// easy peasy, just remove the block with the lowest scores
		Block<S> block = packedBlocks.pollFirst();
		store.freeBlock(block.id);

		size -= block.size;
	}

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
}
