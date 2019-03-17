/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.astar.conf;


import edu.duke.cs.osprey.confspace.ConfSearch;

import java.lang.ref.SoftReference;
import java.math.BigInteger;
import java.util.Iterator;
import java.util.LinkedHashSet;
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
public class ConfSearchCache {

	public class Entry implements ConfSearch {

		private final Supplier<ConfSearch> factory;

		private long numConfs = 0;
		private boolean isExhausted = false;
		private ConfSearch strongRef = null;
		private SoftReference<ConfSearch> softRef = null;

		private Entry(Supplier<ConfSearch> factory) {
			this.factory = factory;
			getOrMakeTree();
		}

		private ConfSearch getOrMakeTree() {

			// check the soft ref to see if we still have a tree
			// (it could have been collected by the GC)
			if (softRef != null) {
				ConfSearch tree = softRef.get();
				if (tree != null) {
					markUsed(tree);
					return tree;
				}
			}

			// don't have a tree, make a new one
			ConfSearch tree = factory.get();

			// and put it back to where it was
			for (int i=0; i<numConfs; i++) {
				tree.nextConf();
			}

			// recently-used entries are always protected from garbage collection
			softRef = new SoftReference<>(tree);
			markUsed(tree);

			return tree;
		}

		private void markUsed(ConfSearch tree) {

			// protect from garbage collection by holding a strong reference
			strongRef = tree;

			// if capacity restrictions are turned on, manage recency and GC protections
			if (minCapacity != null) {

				recentEntries.remove(this);
				recentEntries.add(this);

				// if we're over capacity, expose the least recently used trees to garbage collection
				if (recentEntries.size() > minCapacity) {
					Iterator<Entry> iter = recentEntries.iterator();

					// get rid of the strong reference, so we only have the soft reference
					iter.next().strongRef = null;
					iter.remove();
				}
			}
		}

		public void clearRefs() {
			softRef = null;
			strongRef = null;
		}

		public boolean isProtected() {
			return strongRef != null;
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


	public final Integer minCapacity;

	private final LinkedHashSet<Entry> recentEntries = new LinkedHashSet<>();

	public ConfSearchCache(Integer minCapacity) {
		this.minCapacity = minCapacity;
	}

	public Entry make(Supplier<ConfSearch> factory) {
		return new Entry(factory);
	}
}
