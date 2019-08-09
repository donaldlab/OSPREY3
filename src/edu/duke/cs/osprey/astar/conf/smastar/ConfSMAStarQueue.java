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
