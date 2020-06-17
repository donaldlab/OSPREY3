package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class NodeProcessor implements AutoCloseable {

	public final TaskExecutor nodeTasks;
	public final TaskExecutor minimizationTasks;
	public final SeqDB seqdb;
	public final NodeDB nodedb;
	public final StateInfo[] stateInfos;

	private final List<MinimizationThread> minimizationThreads;

	public NodeProcessor(TaskExecutor nodeTasks, TaskExecutor minimizationTasks, SeqDB seqdb, NodeDB nodedb, StateInfo[] stateInfos) {

		this.nodeTasks = nodeTasks;
		this.minimizationTasks = minimizationTasks;
		this.seqdb = seqdb;
		this.nodedb = nodedb;
		this.stateInfos = stateInfos;

		// start the minimization threads
		minimizationThreads = Arrays.stream(stateInfos)
			.map(stateInfo -> new MinimizationThread(stateInfo, this, minimizationTasks))
			.collect(Collectors.toList());
	}

	@Override
	public void close() {

		// clean up the minimization threads
		for (var thread : minimizationThreads) {
			thread.askToStop();
		}
		for (var thread : minimizationThreads) {
			thread.join();
		}
	}

	private void log(String msg, Object ... args) {
		nodedb.member.log(msg, args);
	}

	public boolean processFor(Directions directions, long duration, TimeUnit timeUnit) {

		boolean foundNodes = false;

		long stopNs = System.nanoTime() + timeUnit.toNanos(duration);
		while (System.nanoTime() < stopNs) {

			var nodeInfo = getNode(directions);
			switch (nodeInfo.result) {

				// process the node in a task
				case GotNode -> {
					foundNodes = true;
					nodeTasks.submit(
						() -> process(nodeInfo),
						ignored -> {}
					);
				}

				// otherwise, wait a bit and try again
				default -> ThreadTools.sleep(100, TimeUnit.MILLISECONDS);
			}
		}
		nodeTasks.waitForFinish();

		return foundNodes;
	}

	public enum Result {

		GotNode(true),
		NoInfo(false),
		NoNode(false);

		public final boolean gotNode;

		Result(boolean gotNode) {
			this.gotNode = gotNode;
		}
	}

	public static class NodeInfo {

		public final Result result;
		public final NodeIndex.Node node;
		public final RCs tree;

		public NodeInfo(Result result, NodeIndex.Node node, RCs tree) {
			this.result = result;
			this.node = node;
			this.tree = tree;
		}

		public NodeInfo(Result result) {
			this(result, null, null);
		}
	}

	public Sequence makeSeq(int statei, int[] conf) {
		if (stateInfos[statei].config.state.isSequenced) {
			return seqdb.confSpace.seqSpace.makeSequence(stateInfos[statei].confSpace, conf);
		} else {
			return null;
		}
	}

	public NodeInfo getNode(Directions directions) {

		// get the currently focused state
		int statei = directions.getFocusedStatei();
		if (statei < 0) {
			return new NodeInfo(Result.NoInfo);
		}

		// get the tree for this state
		RCs tree = directions.getTree(statei);
		if (tree == null) {
			return new NodeInfo(Result.NoInfo);
		}

		// get the next node from that state
		var node = nodedb.removeHigh(statei);
		if (node == null) {
			return new NodeInfo(Result.NoNode);
		}

		// all is well
		return new NodeInfo(Result.GotNode, node, tree);
	}

	public NodeInfo process(NodeInfo nodeInfo) {

		// TODO: NEXTTIME: something in here is racing!!

		var stateInfo = stateInfos[nodeInfo.node.statei];

		var confIndex = stateInfo.makeConfIndex();
		Conf.index(nodeInfo.node.conf, confIndex);

		// is this a leaf node?
		if (confIndex.isFullyDefined()) {

			// yup, add it to the minimization queue
			minimizationThreads.get(nodeInfo.node.statei).addNode(nodeInfo.node);

		} else {

			var seqdbBatch = seqdb.batch();

			// remove the old bound from SeqDB
			seqdbBatch.subZSumUpper(
				stateInfo.config.state,
				makeSeq(nodeInfo.node.statei, nodeInfo.node.conf),
				nodeInfo.node.score
			);

			// pick the next position to expand
			// TODO: allow permuting the positions
			int posi = confIndex.numDefined;

			// expand the node
			for (int confi : nodeInfo.tree.get(posi)) {
				confIndex.assignInPlace(posi, confi);

				// compute an upper bound for the assignment
				var zSumUpper = stateInfo.zSumUpper(confIndex, nodeInfo.tree).normalize(true);

				// update nodedb
				var childNode = new NodeIndex.Node(nodeInfo.node.statei, Conf.make(confIndex), zSumUpper);
				nodedb.addLocal(childNode);

				// add the new bound
				seqdbBatch.addZSumUpper(
					stateInfo.config.state,
					makeSeq(nodeInfo.node.statei, childNode.conf),
					childNode.score
				);

				confIndex.unassignInPlace(posi);
			}

			seqdbBatch.save();
		}

		return nodeInfo;
	}

	public void handleDrops(Stream<NodeIndex.Node> nodes) {

		// NOTE: this is called on the NodeDB thread!

		var batch = seqdb.batch();
		nodes.forEach(node ->
			batch.drop(
				stateInfos[node.statei].config.state,
				makeSeq(node.statei, node.conf),
				node.score
			)
		);
		batch.save();
	}
}
