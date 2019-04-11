package edu.duke.cs.osprey.dof;


import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.tools.HashCalculator;

import java.util.HashMap;
import java.util.Map;


/**
 * A thread-safe cache for MutAlignment instances, which are relatively expensive to compute
 */
public class MutAlignmentCache {

	private static class Key {

		final ResidueTemplate src;
		final ResidueTemplate dst;

		Key(ResidueTemplate src, ResidueTemplate dst) {
			this.src = src;
			this.dst = dst;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				System.identityHashCode(src),
				System.identityHashCode(dst)
			);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof Key && equals((Key)other);
		}

		public boolean equals(Key other) {
			return this.src == other.src
				&& this.dst == other.dst;
		}
	}

	private final Map<Key,MutAlignment> cache = new HashMap<>();

	public synchronized MutAlignment get(ResidueTemplate src, ResidueTemplate dst) {
		return cache.computeIfAbsent(new Key(src, dst), (key) -> new MutAlignment(src, dst));
	}
}
