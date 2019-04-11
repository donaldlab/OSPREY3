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

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ConfSMAStarNode implements ConfAStarNode {

	public static enum State {

		/** never spawned yet */
		Unspawned,

		/** in the child array */
		Spawned,

		/** not spawned now, never need to spawn again */
		Finished,

		/** not spawned now, might need to spawn again */
		Forgotten
	}

	/* SMA* memory usage:

		memory usage per node (according to VisualVM) is 160 bytes + 24 bytes per child
		so 1 GiB of memory could fit ~5.8 M nodes with 1 child each,
			or ~420 K nodes with 100 children each

		the SMA* queue uses roughly 250 bytes per node also
		only a portion of the nodes are in the queue at any one time though
	*/

	public final ConfSMAStarNode parent;
	public final int index;
	public final int depth;
	public final int pos;
	public final int rc;

	private double gscore = 0.0;
	private double hscore = 0.0;
	private double fscore = Double.NaN;

	private ConfSMAStarNode[] spawnedChildren = null;
	private State[] childStates = null;
	private double[] forgottenScores = null;


	/** make the root node */
	public ConfSMAStarNode() {
		this(null, -1, 0, -1, -1);
	}

	/** make a non-root node */
	private ConfSMAStarNode(ConfSMAStarNode parent, int index, int depth, int pos, int rc) {
		this.parent = parent;
		this.index = index;
		this.depth = depth;
		this.pos = pos;
		this.rc = rc;
	}

	private void allocateChildren(int numChildren) {
		spawnedChildren = new ConfSMAStarNode[numChildren];
		Arrays.fill(spawnedChildren, null);
		childStates = new State[numChildren];
		Arrays.fill(childStates, State.Unspawned);
		forgottenScores = new double[numChildren];
		Arrays.fill(forgottenScores, Double.NaN);
	}

	@Override
	public void getConf(int[] conf) {
		ConfSMAStarNode node = this;
		while (node.depth > 0) {
			conf[node.pos] = node.rc;
			node = node.parent;
		}
	}

	@Override
	public double getScore() {
		return fscore;
	}

	public void setScore(double val) {
		fscore = val;
	}

	public double getScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getScore(), optimizer);
	}

	public void setScore(double val, MathTools.Optimizer optimizer) {
		setScore(Tools.optimizeScore(val, optimizer));
	}

	@Override
	public double getGScore() {
		return gscore;
	}

	@Override
	public void setGScore(double val) {
		gscore = val;
	}

	@Override
	public double getHScore() {
		return hscore;
	}

	@Override
	public void setHScore(double val) {
		hscore = val;
	}

	@Override
	public int getLevel() {
		return depth;
	}

	@Override
	public void index(ConfIndex index) {
		index.numDefined = depth;
		ConfSMAStarNode node = this;
		while (node.depth > 0) {
			int i = node.depth - 1;
			index.definedPos[i] = node.pos;
			index.definedRCs[i] = node.rc;
			node = node.parent;
		}
		index.node = this;
		index.sortDefined();
		index.updateUndefined();
	}

	public int getNextChildIndex(int numChildren) {

		if (spawnedChildren == null) {
			allocateChildren(numChildren);
		}

		// pick children we haven't spawned yet first
		for (int i=0; i<childStates.length; i++) {
			if (childStates[i] == State.Unspawned) {
				return i;
			}
		}

		// otherwise, pick the lowest child we've forgotten about
		int bestIndex = -1;
		double bestForgottenScore = Double.POSITIVE_INFINITY;
		for (int i=0; i<childStates.length; i++) {
			if (childStates[i] == State.Forgotten) {
				if (forgottenScores[i] < bestForgottenScore) {
					bestForgottenScore = forgottenScores[i];
					bestIndex = i;
				}
			}
		}

		if (bestIndex >= 0) {
			return bestIndex;
		}

		throw new Error("No more children to spawn");
	}

	@Override
	public ConfAStarNode assign(int pos, int rc) {
		throw new UnsupportedOperationException("need the child index, call spawnChild instead");
	}

	public ConfSMAStarNode spawnChild(int pos, int rc, int index) {

		ConfSMAStarNode child = new ConfSMAStarNode(this, index, depth + 1, pos, rc);
		spawnedChildren[index] = child;
		childStates[index] = State.Spawned;
		forgottenScores[index] = Double.NaN;

		return child;
	}

	public void forgetChild(ConfSMAStarNode child) {
		assert (spawnedChildren[child.index] == child);
		childStates[child.index] = ConfSMAStarNode.State.Forgotten;
		forgottenScores[child.index] = child.fscore;
		spawnedChildren[child.index] = null;
	}

	public int finishChild(ConfSMAStarNode child, ConfSMAStarQueue q) {

		assert (spawnedChildren[child.index] == child);

		// remove the child from the queue, before changing the fscore
		q.removeOrAssert(child);

		// flag this child as finished and update the fscores
		child.fscore = Double.POSITIVE_INFINITY;
		childStates[child.index] = State.Finished;
		forgottenScores[child.index] = Double.POSITIVE_INFINITY;
		backup(q);

		// remove the child from the tree
		spawnedChildren[child.index] = null;
		int numNodesRemoved = 1;

		// recurse up the tree
		ConfSMAStarNode node = this;
		while (node != null && node.allChildrenFinished()) {

			// remove the node from the queue, if needed
			q.remove(node);

			// remove it from the tree too, if needed
			if (node.parent != null) {
				node.parent.childStates[node.index] = State.Finished;
				node.parent.forgottenScores[node.index] = Double.POSITIVE_INFINITY;
				node.parent.spawnedChildren[node.index] = null;
			}
			numNodesRemoved++;

			node = node.parent;
		}

		return numNodesRemoved;
	}

	/** ie, the "backup" operation described by the SMA* paper */
	public void backup(ConfSMAStarQueue q) {

		ConfSMAStarNode node = this;
		while (node != null && node.haveAllChildScores()) {

			double oldScore = node.fscore;

			// update the fscore to the min of any known successor score
			double newScore = Double.POSITIVE_INFINITY;
			for (ConfSMAStarNode child : node.spawnedChildren) {
				if (child != null) {
					newScore = Math.min(newScore, child.fscore);
				}
			}
			for (double score : node.forgottenScores) {
				if (!Double.isNaN(score)) {
					newScore = Math.min(newScore, score);
				}
			}

			// no change? we're done here
			if (newScore == oldScore) {
				break;
			}

			// update the node score (take the node out of the queue while changing the fscore)
			boolean wasRemoved = q.remove(node);
			node.fscore = newScore;
			if (wasRemoved) {
				q.addOrAssert(node);
			}

			node = node.parent;
		}
	}

	public String getPath() {

		assert (depth == 0 || parent != null) : String.format("%d:%d <- null", pos, rc);

		StringBuilder buf = new StringBuilder();
		ConfSMAStarNode node = this;
		if (node.depth == 0) {
			buf.append("(root)");
		}
		while (node.depth > 0) {
			if (buf.length() > 0) {
				buf.append(" <- ");
			}
			buf.append(node.pos);
			buf.append(":");
			buf.append(node.rc);
			assert (node.parent != null) : buf.toString() + " <- null";
			node = node.parent;
		}
		return buf.toString();
	}

	@Override
	public String toString() {
		return String.format("%s   %9.4f (%9.4f)    children %s",
			getPath(),
			fscore,
			forgottenScores == null ? Double.NaN : Arrays.stream(forgottenScores)
				.filter(score -> !Double.isNaN(score))
				.min()
				.orElse(Double.NaN),
			childStates == null ? "(none)" : IntStream.range(0, childStates.length)
				.mapToObj(i -> String.format("%d:%s:%.4f:%.4f", i, childStates[i], forgottenScores[i], spawnedChildren[i] == null ? Double.NaN : spawnedChildren[i].fscore))
				.collect(Collectors.toList())
		);
	}

	/** ie, were all the children spawned at least once, so we've seen their fscores? */
	private boolean haveAllChildScores() {
		for (State state : childStates) {
			if (state == State.Unspawned) {
				return false;
			}
		}
		return true;
	}

	/** ie, is it possible to spawn any children in the future? */
	public boolean canSpawnChildren() {
		for (State state : childStates) {
			if (state == State.Unspawned || state == State.Forgotten) {
				return true;
			}
		}
		return false;
	}

	private boolean allChildrenFinished() {
		for (State state : childStates) {
			if (state != State.Finished) {
				return false;
			}
		}
		return true;
	}

	public boolean hasSpawnedChildren() {
		if (childStates == null) {
			return false;
		}
		for (State state : childStates) {
			if (state == State.Spawned) {
				return true;
			}
		}
		return false;
	}
}
