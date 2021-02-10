package edu.duke.cs.osprey.coffee.nodedb;

import java.nio.ByteBuffer;
import java.util.*;


/**
 * An index that allows queries for items with high scores,
 * but drops items with low scores to make extra space.
 *
 * However, since multiple indices can share the same memory,
 * it may not always be possible to free up space in this index.
 *
 * Implementation are NOT thread-safe unless otherwise specified.
 */
public interface FixedIndex<S extends Comparable<S>, T extends FixedIndex.Indexable<S>> {

	interface Indexable<S> {
		S score();
	}

	interface Serializer<T> {
		int bytes();
		void serialize(ByteBuffer out, T data);
		T deserialize(ByteBuffer in);
	}

	/** Returns the max number of items per block. */
	int blockCapacity();

	/** Returns the number of items currently stored. */
	long size();

	/** Returns the number of items that can be added without freeing up space. */
	long freeSpace();

	/**
	 * Adds an item to the index.
	 * Returns true if the item was added, false if there wasn't enough space.
	 */
	boolean add(T item);

	/** Returns the highest score of any currently stored item. */
	S highestScore();

	/** Removes the highest-scoring item currently stored. */
	T removeHighest();

	/** Drops low-scoring items to make more free space, and adds them to the dropped queue. */
	void freeUpSpace();

	/** Returns the queue of dropped items. */
	Deque<T> dropped();

	/** Removes all items. */
	void clear();
}
