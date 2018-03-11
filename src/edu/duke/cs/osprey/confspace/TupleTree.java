package edu.duke.cs.osprey.confspace;


import edu.duke.cs.osprey.tools.UnpossibleError;

import java.util.ArrayList;
import java.util.List;
import java.util.function.BiConsumer;


/**
 * efficient storage for high-dimensional RC tuples based on a prefix tree (aka "trie")
 * @param <T> data type to store at the tuple
 */
public class TupleTree<T> {

	private static final int NotATuple = -1;

	private class Node {

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

		@SuppressWarnings("unchecked")
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
				if (baseTuple.pos.get(i) != tuple.pos.get(i) || baseTuple.RCs.get(i) != tuple.RCs.get(i)) {
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

	public void forEachIn(int[] conf, BiConsumer<RCTuple,T> callback) {

		// check the base tuple if needed
		if (baseTuple != null) {
			for (int i=0; i<baseTuple.size(); i++) {

				int pos = baseTuple.pos.get(i);
				int rc = baseTuple.RCs.get(i);

				if (conf[pos] != rc) {
					return;
				}
			}
		}

		forEachIn(conf, root, conf.length, callback);
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
}
