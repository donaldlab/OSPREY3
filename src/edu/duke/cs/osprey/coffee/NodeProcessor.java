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
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.Stopwatch;

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
	private final List<int[]> posPermutations;

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

		// pick an ordering of the positions that speeds things up
		// for now, sequence positions at the top is a good heuristic
		// TODO: within seq/conf strata, high branch factor at the top?
		posPermutations = Arrays.stream(stateInfos)
			.map(stateInfo -> Arrays.stream(stateInfo.confSpace.positions)
				.sorted(Comparator.comparing(pos -> !pos.hasMutations))
				.mapToInt(pos -> pos.index)
				.toArray()
			)
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

		var statei = nodeInfo.node.statei;
		var stateInfo = stateInfos[statei];

		var confIndex = stateInfo.makeConfIndex();
		Conf.index(nodeInfo.node.conf, confIndex);

		// is this a leaf node?
		if (confIndex.isFullyDefined()) {

			// yup, see if this node's score is roughly as good as current predictions
			var currentScore = nodedb.perf.score(statei, nodeInfo.node.conf, nodeInfo.node.zSumUpper);
			if (nodeInfo.node.score.exp - currentScore.exp > 1) {

				// TEMP
				int diff = nodeInfo.node.score.exp - currentScore.exp;

				// nope, it should have a much worse score
				// re-score it and put it back into nodedb
				nodedb.add(new NodeIndex.Node(statei, nodeInfo.node.conf, nodeInfo.node.zSumUpper, currentScore));

				/* TEMP
				log("rescore node: %s  %s -> %s   diff %d %b",
					Arrays.toString(nodeInfo.node.conf),
					nodeInfo.node.score,
					currentScore,
					diff, Math.abs(currentScore.exp - nodeInfo.node.score.exp) > 1
				);
				*/

			} else {

				// the score look good, add it to the minimization queue
				minimizationThreads.get(statei).addNode(nodeInfo.node);

			}
		} else {

			// get timing info for the node expansion
			Stopwatch stopwatch = new Stopwatch().start();

			// remove the old bound from SeqDB
			var seqdbBatch = seqdb.batch();
			seqdbBatch.subZSumUpper(
				stateInfo.config.state,
				makeSeq(statei, nodeInfo.node.conf),
				nodeInfo.node.zSumUpper
			);

			// track the reduction in uncertainty for this node
			var reduction = new BigExp(nodeInfo.node.zSumUpper);

			// pick the next position to expand, according to the position permutation
			int posi = posPermutations.get(statei)[confIndex.numDefined];

			// expand the node at the picked position
			for (int confi : nodeInfo.tree.get(posi)) {
				confIndex.assignInPlace(posi, confi);

				// compute an upper bound for the assignment
				var zSumUpper = stateInfo.zSumUpper(confIndex, nodeInfo.tree).normalize(true);

				// update nodedb
				var conf = Conf.make(confIndex);
				var childNode = new NodeIndex.Node(
					statei, conf, zSumUpper,
					nodedb.perf.score(statei, conf, zSumUpper)
				);
				nodedb.add(childNode);

				// add the new bound
				seqdbBatch.addZSumUpper(
					stateInfo.config.state,
					makeSeq(statei, conf),
					zSumUpper
				);

				// update the reduction calculation
				reduction.sub(zSumUpper);

				confIndex.unassignInPlace(posi);
			}

			// compute the uncertainty reduction and update the node performance
			stopwatch.stop();
			nodedb.perf.update(nodeInfo.node, stopwatch.getTimeNs(), reduction);

			/* TEMP
			log("node expansion:   r %s  %10s  r/s %s",
				reduction,
				stopwatch.getTime(2),
				Log.formatBigEngineering(seqdb.bigMath()
					.set(reduction)
					.div(stopwatch.getTimeNs())
					.get()
				)
			);
			*/

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
				node.zSumUpper
			)
		);
		batch.save();
	}
}
