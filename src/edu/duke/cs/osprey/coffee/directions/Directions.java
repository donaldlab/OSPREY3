package edu.duke.cs.osprey.coffee.directions;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.ClusterMember;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;


public class Directions {

	static final String ServiceName = "Coffee";

	public final MultiStateConfSpace confSpace;
	public final ClusterMember member;

	private final AtomicBoolean isRunning = new AtomicBoolean(true);
	private final AtomicInteger focusedStatei = new AtomicInteger(-1);
	private final AtomicReference<RCs[]> trees = new AtomicReference<>(null);
	private final List<Set<Sequence>> finishedSeqs;

	public Directions(MultiStateConfSpace confSpace, ClusterMember member) {

		this.confSpace = confSpace;
		this.member = member;

		finishedSeqs = confSpace.sequencedStates.stream()
			.map(state -> new HashSet<Sequence>())
			.collect(Collectors.toList());

		// register with hazelcast
		member.registerService(ServiceName, this);
	}

	/**
	 * Tell all cluster members to stop processing.
	 */
	public void stop() {
		receiveStop();
		member.sendToOthers(() -> new StopOperation());
	}

	void receiveStop() {
		isRunning.set(false);
	}

	public boolean isRunning() {
		return isRunning.get();
	}

	/**
	 * Tell all cluster members to process nodes in the specified state.
	 */
	public void focus(int statei) {
		receiveFocus(statei);
		member.sendToOthers(() -> new FocusOperation(statei));
	}

	void receiveFocus(int statei) {
		focusedStatei.set(statei);
	}

	public int getFocusedStatei() {
		return focusedStatei.get();
	}

	/**
	 * Tell all cluster members the node tree topology for each state.
	 */
	public void setTrees(RCs[] trees) {
		receiveTrees(trees);
		member.sendToOthers(() -> new TreesOperation(trees));
	}

	void receiveTrees(RCs[] trees) {
		this.trees.set(trees);
	}

	public RCs getTree(int statei) {
		return trees.get()[statei];
	}

	public RCs getTreeOrThrow(int statei) {
		RCs tree = getTree(statei);
		if (tree != null) {
			return tree;
		}
		throw new IllegalStateException("no node tree set for state " + statei);
	}

	/**
	 * Adds the sequences to the finished set for the given state.
	 * newlyFinishedSeqs will be cleared.
	 * Any sequences that were not already finished will be added to newlyFinishedSeqs.
	 */
	public void finishSequences(int sequencedStatei, Set<Sequence> finishedSeqs, Set<Sequence> newlyFinishedSeqs) {

		newlyFinishedSeqs.clear();
		newlyFinishedSeqs.addAll(finishedSeqs);

		var finished = this.finishedSeqs.get(sequencedStatei);
		synchronized (finished) {
			newlyFinishedSeqs.removeAll(finished);
			finished.addAll(newlyFinishedSeqs);
		}

		// if there are any newly added sequences, send those to the other nodes
		if (!newlyFinishedSeqs.isEmpty()) {
			member.sendToOthers(() -> new FinishedSequencesOperation(sequencedStatei, newlyFinishedSeqs));
		}
	}

	void receiveFinishedSequences(int sequencedStatei, int[][] seqs) {
		var finished = finishedSeqs.get(sequencedStatei);
		synchronized (finished) {
			for (int[] seq : seqs) {
				finished.add(new Sequence(confSpace.seqSpace, seq));
			}
		}
	}

	public boolean isFinished(int sequencedStatei, Sequence seq) {
		var finished = finishedSeqs.get(sequencedStatei);
		synchronized (finished) {
			return finished.contains(seq);
		}
	}
}
