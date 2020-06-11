package edu.duke.cs.osprey.coffee.directions;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.ClusterMember;


public class Directions {

	static final String ServiceName = "Coffee";

	public final ClusterMember member;

	private boolean isRunning = true;
	private int focusedStatei = -1;
	private RCs[] trees = null;

	public Directions(ClusterMember member) {

		this.member = member;

		// register with hazelcast
		member.registerService(ServiceName, this);
	}

	/**
	 * Tell all cluster members to stop processing.
	 */
	public void stop() {
		member.sendToOthers(() -> new StopOperation());
	}

	void receiveStop() {
		isRunning = false;
	}

	public boolean isRunning() {
		return isRunning;
	}

	/**
	 * Tell all cluster members to process nodes in the specified state.
	 */
	public void focus(int statei) {
		receiveFocus(statei);
		member.sendToOthers(() -> new FocusOperation(statei));
	}

	void receiveFocus(int statei) {
		focusedStatei = statei;
	}

	public int getFocusedStatei() {
		return focusedStatei;
	}

	/**
	 * Tell all cluster members the node tree topology for each state.
	 */
	public void setTrees(RCs[] trees) {
		receiveTrees(trees);
		member.sendToOthers(() -> new TreesOperation(trees));
	}

	void receiveTrees(RCs[] trees) {
		this.trees = trees;
	}

	public RCs getTree(int statei) {
		return trees[statei];
	}

	public RCs getTreeOrThrow(int statei) {
		RCs tree = getTree(statei);
		if (tree != null) {
			return tree;
		}
		throw new IllegalStateException("no node tree set for state " + statei);
	}

	// TODO: ignore sequences
}
