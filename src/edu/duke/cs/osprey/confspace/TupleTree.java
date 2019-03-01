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

package edu.duke.cs.osprey.confspace;


import edu.duke.cs.osprey.tools.UnpossibleError;

import java.io.Serializable;
import java.util.*;
import java.util.function.BiConsumer;


/**
 * efficient storage for high-dimensional RC tuples based on a prefix tree (aka "trie")
 * @param <T> data type to store at the tuple
 */
public class TupleTree<T> implements Serializable {

	private static final long serialVersionUID = -8048566743262701431L;

	private class Node implements Serializable {

		private static final long serialVersionUID = -47314393953695943L;

		final int pos;
		final int rc;

		/*
			use a list of lists here to act as a lookup table for child nodes
			this is stupidly fast, but can potentially use a lot of memory
			if we need more efficient space usage, we could use hash tables instead of arrays
			but I tried hash tables already and they were slower by a noticeable amount
		 */
		List<List<Node>> children = null;

		RCTuple tuple = null;
		T data = null;

		Node(int pos, int rc) {
			this.pos = pos;
			this.rc = rc;
		}

		public Node makeChild(int pos, int rc) {

			// make sure children exists
			if (children == null) {
				children = new ArrayList<>();
			}

			// make sure children[pos] exists
			while (children.size() <= pos) {
				children.add(null);
			}
			List<Node> nodes = children.get(pos);
			if (nodes == null) {
				nodes = new ArrayList<>();
				children.set(pos, nodes);
			}

			// make sure children[pos][rc] exists
			while (nodes.size() <= rc) {
				nodes.add(null);
			}
			Node child = new Node(pos, rc);
			nodes.set(rc, child);

			return child;
		}

		public Node get(int pos, int rc) {
			if (children == null) {
				return null;
			}
			if (pos >= children.size()) {
				return null;
			}
			List<Node> nodes = children.get(pos);
			if (nodes == null) {
				return null;
			}
			if (rc >= nodes.size()) {
				return null;
			}
			return nodes.get(rc);
		}
	}

	/**
	 * a common tuple prefix shared by all tuples in this tree
	 * useful for implementing a pairwise matrix of tuple trees, for e.g. AbstractTupleMatrix
	 */
	public final RCTuple baseTuple;

	public final Node root = new Node(-1, -1);

	public TupleTree() {
		this(null);
	}

	public TupleTree(RCTuple baseTuple) {

		if (baseTuple != null) {
			baseTuple.checkSortedPositions();
		}

		this.baseTuple = baseTuple;
	}

	private void checkTuple(RCTuple tuple) {

		// make sure positions are sorted
		tuple.checkSortedPositions();

		if (baseTuple != null) {

			// make sure it matches the base
			for (int i=0; i<baseTuple.size(); i++) {
				int bpos = baseTuple.pos.get(i);
				int brc = baseTuple.RCs.get(i);
				int tpos = tuple.pos.get(i);
				int trc = tuple.RCs.get(i);
				if (bpos != tpos || brc != trc) {
					throw new IllegalArgumentException("Tuple " + tuple + " doesn't match base tuple " + baseTuple + " for this tree");
				}
			}
		}
	}

	private int getFirstIndex() {
		if (baseTuple == null) {
			return 0;
		}
		return baseTuple.size();
	}

	public T get(RCTuple tuple) {

		// just in case...
		checkTuple(tuple);

		Node parent = root;
		int firstIndex = getFirstIndex();
		for (int i=tuple.size()-1; i>=firstIndex; i--) {

			int pos = tuple.pos.get(i);
			int rc = tuple.RCs.get(i);
			Node node = parent.get(pos, rc);
			if (node == null) {
				return null;
			}

			if (i == firstIndex) {
				return node.data;
			} else {
				parent = node;
			}
		}

		throw new UnpossibleError();
	}

	public void put(RCTuple tuple, T data) {

		// just in case...
		checkTuple(tuple);

		// hang tuples by their toes!
		// i.e., by greatest pos first
		// conf pos iteration usually follows pos3 < pos2 < pos1 etc
		// pos1, the largest value, is always in the outer loop and should correspond to the top node
		Node parent = root;
		int firstIndex = getFirstIndex();
		for (int i=tuple.size()-1; i>=firstIndex; i--) {

			int pos = tuple.pos.get(i);
			int rc = tuple.RCs.get(i);
			Node node = parent.get(pos, rc);
			if (node == null) {
				node = parent.makeChild(pos, rc);
			}

			if (i == firstIndex) {
				node.data = data;
				node.tuple = tuple;
			} else {
				parent = node;
			}
		}
	}

	public void clear() {
		root.children = null;
		root.data = null;
	}

	private boolean matchesBaseTuple(int[] conf) {

		if (baseTuple == null) {
			return true;
		}

		for (int i=0; i<baseTuple.size(); i++) {

			int pos = baseTuple.pos.get(i);
			int rc = baseTuple.RCs.get(i);

			if (rc != conf[pos]) {
				return false;
			}
		}

		return true;
	}

	public void forEachIn(int[] conf, BiConsumer<RCTuple,T> callback) {
		if (matchesBaseTuple(conf)) {
			forEachIn(conf, root, conf.length, callback);
		}
	}

	private void forEachIn(int[] conf, Node parent, int untilPos, BiConsumer<RCTuple,T> callback) {

		for (int pos=0; pos<untilPos; pos++) {

			// skip unassigned positions
			int rc = conf[pos];
			if (rc == Conf.Unassigned) {
				continue;
			}

			// get the node for this rc, if any
			Node node = parent.get(pos, rc);
			if (node == null) {
				continue;
			}

			// callback if there's a tuple at this node
			if (node.tuple != null) {
				callback.accept(node.tuple, node.data);
			}

			// recurse if possible
			if (node.children != null) {
				forEachIn(conf, node, pos, callback);
			}
		}
	}

	public void forEachIn(int[] conf, int pos1, BiConsumer<RCTuple,T> callback) {
		if (matchesBaseTuple(conf)) {
			forEachIn(conf, pos1, root, pos1 + 1, callback);
		}
	}

	private void forEachIn(int[] conf, int pos1, Node parent, int untilPos, BiConsumer<RCTuple,T> callback) {

		for (int pos=pos1; pos<untilPos; pos++) {

			// skip unassigned positions
			int rc = conf[pos];
			if (rc == Conf.Unassigned) {
				continue;
			}

			// get the node for this rc, if any
			Node node = parent.get(pos, rc);
			if (node == null) {
				continue;
			}

			// callback if there's a tuple at this node that matches pos1
			if (node.tuple != null && node.tuple.pos.contains(pos1)) {
				callback.accept(node.tuple, node.data);
			}

			// recurse if possible
			if (node.children != null) {
				forEachIn(conf, node, pos, callback);
			}
		}
	}

	public void forEachIn(int[] conf, int pos1, int pos2, BiConsumer<RCTuple,T> callback) {
		if (matchesBaseTuple(conf)) {
			forEachIn(conf, pos1, pos2, root, pos1 + 1, callback);
		}
	}

	private void forEachIn(int[] conf, int pos1, int pos2, Node parent, int untilPos, BiConsumer<RCTuple,T> callback) {

		for (int pos=Math.max(pos1, pos2); pos<untilPos; pos++) {

			// skip unassigned positions
			int rc = conf[pos];
			if (rc == Conf.Unassigned) {
				continue;
			}

			// get the node for this rc, if any
			Node node = parent.get(pos, rc);
			if (node == null) {
				continue;
			}

			// callback if there's a tuple at this node that matches pos1, pos2
			if (node.tuple != null && node.tuple.pos.contains(pos1) && node.tuple.pos.contains(pos2)) {
				callback.accept(node.tuple, node.data);
			}

			// recurse if possible
			if (node.children != null) {
				forEachIn(conf, node, pos, callback);
			}
		}
	}

	// TODO: make an iterator instead of writing to a list?
	// TODO: make a forEach style callback?

	public List<RCTuple> makeTuplesList() {
		List<RCTuple> tuples = new ArrayList<>();
		addTuplesToList(tuples, root);
		return tuples;
	}

	private void addTuplesToList(List<RCTuple> tuples, Node node) {

		// any tuples here?
		if (node.tuple != null) {
			tuples.add(node.tuple);
		}

		// recurse
		if (node.children == null) {
			return;
		}
		for (List<Node> nodes : node.children) {
			if (nodes == null) {
				continue;
			}
			for (Node child : nodes) {
				if (child == null) {
					continue;
				}
				addTuplesToList(tuples, child);
			}
		}
	}
}
