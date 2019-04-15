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
