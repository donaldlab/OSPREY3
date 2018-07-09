package edu.duke.cs.osprey.astar.conf.smastar;

import java.util.*;


public class ConfSMAStarQueue {

	private static class Nodes {

		final double score;
		final TreeMap<Integer,ConfSMAStarNode> nodesByDepth = new TreeMap<>();

		Nodes(double score) {
			this.score = score;
		}
	}

	private final TreeMap<Double,Nodes> nodesByScore = new TreeMap<>();


	public boolean add(ConfSMAStarNode node) {

		assert (Double.isFinite(node.getScore()));

		Nodes nodes = getOrMakeNodes(node.getScore());
		if (nodes.nodesByDepth.containsKey(node.depth)) {
			return false;
		}
		nodes.nodesByDepth.put(node.depth, node);
		return true;
	}

	public void addOrAssert(ConfSMAStarNode node) {
		boolean wasAdded = add(node);
		assert (wasAdded);
	}

	private Nodes getOrMakeNodes(double score) {
		Nodes nodes = nodesByScore.get(score);
		if (nodes == null) {
			nodes = new Nodes(score);
			nodesByScore.put(score, nodes);
		}
		return nodes;
	}

	public boolean remove(ConfSMAStarNode node) {

		Nodes nodes = nodesByScore.get(node.getScore());
		if (nodes == null) {
			return false;
		}

		if (nodes.nodesByDepth.get(node.depth) != node) {
			return false;
		}

		ConfSMAStarNode removedNode = nodes.nodesByDepth.remove(node.depth);
		assert (removedNode == node);

		// remove the nodes too, if empty
		if (nodes.nodesByDepth.isEmpty()) {
			nodesByScore.remove(nodes.score);
		}

		return true;
	}

	public void removeOrAssert(ConfSMAStarNode node) {
		boolean wasRemoved = remove(node);
		assert (wasRemoved);
	}

	public boolean isEmpty() {
		return nodesByScore.isEmpty();
	}

	public ConfSMAStarNode getLowestDeepest() {

		if (isEmpty()) {
			return null;
		}

		Nodes lowestNodes = nodesByScore.firstEntry().getValue();
		return lowestNodes.nodesByDepth.lastEntry().getValue();
	}

	public ConfSMAStarNode removeHighestShallowestLeaf() {

		if (isEmpty()) {
			throw new NoSuchElementException("queue is empty");
		}

		Nodes highestNodes = nodesByScore.lastEntry().getValue();
		ConfSMAStarNode leaf = null;
		Iterator<ConfSMAStarNode> iter = highestNodes.nodesByDepth.values().iterator();
		while (iter.hasNext()) {
			ConfSMAStarNode node = iter.next();
			if (!node.hasSpawnedChildren()) {
				leaf = node;
				iter.remove();
				break;
			}
		}

		// remove the nodes, if needed
		if (highestNodes.nodesByDepth.isEmpty()) {
			nodesByScore.remove(highestNodes.score);
		}

		return leaf;
	}
}
