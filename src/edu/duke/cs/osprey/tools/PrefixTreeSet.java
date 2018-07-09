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

package edu.duke.cs.osprey.tools;

import java.util.*;

public class PrefixTreeSet<T> {
	// TODO: implement set interface?

	// probably not the most efficient implementation, but it'll do for now

	private class Node {

		public final T value;
		public final boolean isLeaf;
		public Map<T,Node> children = new HashMap<>(); // TODO: can optimize by making null by default

		public Node(T value, boolean isLeaf) {
			this.value = value;
			this.isLeaf = isLeaf;
		}
	}

	private Map<T,Node> roots = new HashMap<>();

	public void add(List<T> sequence) {

		Map<T,Node> nodes = roots;

		Iterator<T> iter = sequence.iterator();
		while (iter.hasNext()) {
			T item = iter.next();
			boolean isLastItem = !iter.hasNext();

			Node node = nodes.get(item);
			if (node == null) {
				node = new Node(item, isLastItem);
				nodes.put(item, node);
			}

			nodes = node.children;
		}
	}

	public boolean contains(List<T> sequence) {

		if (sequence.isEmpty()) {
			return false;
		}

		Map<T,Node> nodes = roots;

		Iterator<T> iter = sequence.iterator();
		while (iter.hasNext()) {
			T item = iter.next();
			boolean isLastItem = !iter.hasNext();

			Node node = nodes.get(item);
			if (node == null) {
				return false;
			}

			if (isLastItem) {
				return node.isLeaf;
			}

			nodes = node.children;
		}

		// execution should never reach this line, but the compiler can't tell that
		throw new Error("unpossible");
	}
}
