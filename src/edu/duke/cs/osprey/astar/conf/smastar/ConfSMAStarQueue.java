package edu.duke.cs.osprey.astar.conf.smastar;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

public class ConfSMAStarQueue implements Iterable<ConfSMAStarNode> {

	private List<ConfSMAStarNode> nodes = new ArrayList<>();

	public void add(ConfSMAStarNode node) {

		assert (Double.isFinite(node.fscore));

		nodes.add(node);
	}

	public void remove(ConfSMAStarNode node) {

		boolean wasRemoved = nodes.remove(node);
		assert (wasRemoved);
	}

	public boolean contains(ConfSMAStarNode node) {
		return nodes.contains(node);
	}

	@Override
	public Iterator<ConfSMAStarNode> iterator() {
		return nodes.iterator();
	}

	public ConfSMAStarNode getLowestDeepest() {

		assert (!nodes.isEmpty());

		// first, find the lowest nodes
		double lowestScore = nodes.stream()
			.mapToDouble(node -> node.fscore)
			.min()
			.getAsDouble();
		List<ConfSMAStarNode> lowestNodes = nodes.stream()
			.filter(node -> node.fscore == lowestScore)
			.collect(Collectors.toList());
		assert (!lowestNodes.isEmpty());

		// then find the deepest node
		int deepestDepth = lowestNodes.stream()
			.mapToInt(node -> node.depth)
			.max()
			.getAsInt();
		List<ConfSMAStarNode> lowestDeepestNodes = lowestNodes.stream()
			.filter(node -> node.depth == deepestDepth)
			.collect(Collectors.toList());
		return lowestDeepestNodes.get(0);
	}

	public ConfSMAStarNode removeHighestShallowestLeaf() {

		assert (!nodes.isEmpty());

		// first, find the highest nodes
		double highestScore = nodes.stream()
			.mapToDouble(node -> node.fscore)
			.max()
			.getAsDouble();
		List<ConfSMAStarNode> highestNodes = nodes.stream()
			.filter(node -> node.fscore == highestScore)
			.collect(Collectors.toList());
		assert (!highestNodes.isEmpty());

		// then, find the leaf nodes
		List<ConfSMAStarNode> leafNodes = highestNodes.stream()
			.filter(node -> !node.hasSpawnedChildren())
			.collect(Collectors.toList());
		assert (!leafNodes.isEmpty()) : "no leaves\n\t" + String.join("\n\t", highestNodes.stream().map(node -> node.toString()).collect(Collectors.toList()));

		// then, find the shallowest node
		int shallowestDepth = leafNodes.stream()
			.mapToInt(node -> node.depth)
			.min()
			.getAsInt();
		List<ConfSMAStarNode> highestShallowestNodes = highestNodes.stream()
			.filter(node -> node.depth == shallowestDepth)
			.collect(Collectors.toList());
		ConfSMAStarNode node = highestShallowestNodes.get(0);

		remove(node);

		return node;
	}

	public int size() {
		return nodes.size();
	}

	public boolean isEmpty() {
		return nodes.isEmpty();
	}

	public String dumpScores() {
		return nodes.stream()
			.sorted(Comparator.comparing(node -> node.fscore))
			.map(node -> String.format("%.4f", node.fscore))
			.collect(Collectors.toList())
			.toString();
	}
}
