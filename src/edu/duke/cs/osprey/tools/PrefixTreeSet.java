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
