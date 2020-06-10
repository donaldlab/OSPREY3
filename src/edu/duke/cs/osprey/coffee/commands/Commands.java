package edu.duke.cs.osprey.coffee.commands;

import edu.duke.cs.osprey.coffee.ClusterMember;


public class Commands {

	static final String ServiceName = "Coffee";

	public final ClusterMember member;

	private boolean isRunning = true;
	private int focusedStatei = -1;

	public Commands(ClusterMember member) {

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
		focusedStatei = statei;
		member.sendToOthers(() -> new FocusOperation(statei));
	}

	void receiveFocus(int statei) {
		focusedStatei = statei;
	}

	public int getFocusedStatei() {
		return focusedStatei;
	}

	// TODO: ignore sequences
}
