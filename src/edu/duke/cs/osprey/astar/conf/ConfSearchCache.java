package edu.duke.cs.osprey.astar.conf;


import edu.duke.cs.osprey.confspace.ConfSearch;

import java.lang.ref.SoftReference;
import java.math.BigInteger;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Supplier;


/**
 * A cache for ConfSearch instances where the most N recently used
 * instances are protected by strong references and cannot be garbage
 * collected. The remaining instances are held by soft references and
 * will be garbage collected when running low on heap space.
 *
 * Collected trees will be re-instantiated and enumerated to their
 * last known position when accessed again.
 */
public class ConfSearchCache<T> {

	private static class Entry implements ConfSearch {

		final Supplier<ConfSearch> factory;

		long numConfs = 0;
		boolean isExhausted = false;
		ConfSearch strongRef = null;
		SoftReference<ConfSearch> softRef = null;

		public Entry(Supplier<ConfSearch> factory) {
			this.factory = factory;
		}

		ConfSearch getOrMakeTree() {

			// check the soft ref to see if we still have a tree
			// (it could have been collected by the GC)
			if (softRef != null) {
				ConfSearch tree = softRef.get();
				if (tree != null) {
					return tree;
				}
			}

			// don't have a tree, make a new one
			ConfSearch tree = factory.get();

			// and put it back to where it was
			for (int i=0; i<numConfs; i++) {
				tree.nextConf();
			}

			// keep only a soft reference by default
			softRef = new SoftReference<>(tree);

			return tree;
		}

		void protect() {
			strongRef = getOrMakeTree();
		}

		void expose() {
			strongRef = null;
		}

		void clearRefs() {
			softRef = null;
			strongRef = null;
		}

		@Override
		public BigInteger getNumConformations() {
			return getOrMakeTree().getNumConformations();
		}

		@Override
		public ScoredConf nextConf() {

			// no more confs? don't bother with the tree
			if (isExhausted) {
				return null;
			}

			// get the next conf
			ScoredConf conf = getOrMakeTree().nextConf();

			// and keep track of which conf we're on
			if (conf == null) {
				isExhausted = true;

				// and let GC take the tree
				clearRefs();

			} else {
				numConfs++;
			}

			return conf;
		}
	}


	public final int minCapacity;

	private final Map<T,Entry> entries = new HashMap<>();
	private final Deque<Entry> recentEntries = new ArrayDeque<>();

	public ConfSearchCache(int minCapacity) {
		this.minCapacity = minCapacity;
	}

	public ConfSearch getOrMake(T key, Supplier<ConfSearch> factory) {

		// get the entry, or make a new one if needed
		Entry entry = entries.get(key);
		if (entry == null) {
			entry = new Entry(factory);
			entries.put(key, entry);
		}

		// put this entry in the protected capacity
		entry.protect();
		if (!recentEntries.contains(entry)) {
			recentEntries.offer(entry);
		}

		// if we're over capacity, expose the extra trees to garbage collection
		while (recentEntries.size() > minCapacity) {
			recentEntries.pop().expose();
		}

		return entry;
	}

	/**
	 * Clears all references to ConfSearch instances,
	 * forcing them to be re-instantiated on next access.
	 *
	 * Probably only useful for testing (ie simulating garbage collection under memory pressure).
	 */
	public void clearRefs() {
		for (Entry entry : entries.values()) {
			entry.clearRefs();
		}
	}
}
