package edu.duke.cs.osprey.astar.conf.smastar;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ConfSMAStarNode implements ConfAStarNode {

	public static enum State {
		Unspawned, // never spawned yet
		Spawned, // in the child array
		Finished, // not spawned now, never need to spawn again
		Forgotten // not spawned now, might need to spawn again
	}

	// TODO: encapsulation, finality

	public ConfSMAStarNode parent = null;
	public int index = -1;
	public ConfSMAStarNode[] spawnedChildren = null;
	public State[] childStates = null;
	public double[] forgottenScores = null;
	public int depth = 0;
	public double gscore = 0.0;
	public double hscore = 0.0;
	public double fscore = Double.NaN;
	public int pos = -1;
	public int rc = -1;

	public void allocateChildren(int numChildren) {
		spawnedChildren = new ConfSMAStarNode[numChildren];
		Arrays.fill(spawnedChildren, null);
		childStates = new State[numChildren];
		Arrays.fill(childStates, State.Unspawned);
		forgottenScores = new double[numChildren];
		Arrays.fill(forgottenScores, Double.NaN);
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

	public ConfSMAStarNode spawnChild(int pos, int rc, int index, double gscore, double hscore) {
		ConfSMAStarNode child = new ConfSMAStarNode();
		child.parent = this;
		child.index = index;
		child.depth = depth + 1;
		child.gscore = gscore;
		child.hscore = hscore;
		child.fscore = gscore + hscore;
		child.pos = pos;
		child.rc = rc;
		spawnedChildren[index] = child;

		// reset the forgotten score if needed
		if (childStates[index] == State.Forgotten) {
			forgottenScores[index] = Double.NaN;
		}

		childStates[index] = State.Spawned;
		return child;
	}

	public boolean hasSpawnableChildren() {
		for (State state : childStates) {
			if (state == State.Unspawned || state == State.Forgotten) {
				return true;
			}
		}
		return false;
	}

	public int getNextIndex(int numChildren) {

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

	public void removeFromParent() {

		assert (parent.spawnedChildren[index] == this);

		parent.spawnedChildren[index] = null;
		this.parent = null;
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
			assert (node.depth == 0 || node.parent != null) : buf.toString() + " <- null";
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

	@Override
	public ConfAStarNode assign(int pos, int rc) {
		throw new UnsupportedOperationException();
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

	@Override
	public double getGScore() {
		return gscore;
	}

	@Override
	public void setGScore(double val) {
		throw new UnsupportedOperationException();
	}

	@Override
	public double getHScore() {
		return hscore;
	}

	@Override
	public void setHScore(double val) {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getLevel() {
		throw new UnsupportedOperationException();
	}

	public boolean hasUnspawnedChildren() {
		for (State state : childStates) {
			if (state == State.Unspawned) {
				return true;
			}
		}
		return false;
	}

	public void backup(ConfSMAStarQueue q, int maxDepth) {

		if (hasUnspawnedChildren()) {
			return;
		}

		double oldScore = fscore;

		// update the fscore to the min of any known successor score
		double newScore = Double.POSITIVE_INFINITY;
		for (ConfSMAStarNode child : spawnedChildren) {
			if (child != null) {
				newScore = Math.min(newScore, child.fscore);
			}
		}
		for (double score : forgottenScores) {
			if (!Double.isNaN(score)) {
				newScore = Math.min(newScore, score);
			}
		}

		if (newScore != oldScore) {

			// TODO: re-index after score change
			fscore = newScore;

			// recurse if possible
			if (parent != null) {
				parent.backup(q, maxDepth);
			}
		}
	}

	public double getMinForgottenScore() {
		if (forgottenScores == null) {
			return Double.NaN;
		}
		double minScore = Double.NaN;
		for (double score : forgottenScores) {
			if (!Double.isNaN(score)) {
				if (Double.isNaN(minScore) || score < minScore) {
					minScore = score;
				}
			}
		}
		return minScore;
	}

	public boolean allChildrenSpawned() {
		for (State state : childStates) {
			if (state == State.Unspawned || state == State.Forgotten) {
				return false;
			}
		}
		return true;
	}

	public boolean allChildrenFinished() {
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
