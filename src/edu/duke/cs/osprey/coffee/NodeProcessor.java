package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.SeqDB;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class NodeProcessor implements AutoCloseable {

	private static class MinimizationQueue {

		final int size;
		final Deque<NodeIndex.Node> queue;

		MinimizationQueue(StateInfo stateInfo, int size) {
			this.size = size;
			queue = new ArrayDeque<>(size);
		}

		/**
		 * Adds the node to the queue.
		 * If the queue is full, drains the queue and returns the nodes.
		 */
		public synchronized List<NodeIndex.Node> add(NodeIndex.Node node) {

			assert (queue.size() < size);
			queue.add(node);

			if (queue.size() == size) {
				var nodes = new ArrayList<>(queue);
				queue.clear();
				return nodes;
			}

			return null;
		}

		/**
		 * Removes all nodes from the queue.
		 */
		public synchronized List<NodeIndex.Node> flush() {

			if (!queue.isEmpty()) {
				var nodes = new ArrayList<>(queue);
				queue.clear();
				return nodes;
			}

			return null;
		}
	}

	private class NodeThread extends Thread {

		final int id;
		final Directions directions;

		NodeThread(int id, Directions directions) {

			this.id = id;
			this.directions = directions;

			setName("Node-" + id);
			setDaemon(true);
			start();
		}

		@Override
		public void run() {
			while (directions.isRunning()) {

				// get a node
				var nodeInfo = getNode(directions);
				if (!nodeInfo.result.gotNode) {
					// no nodes, wait a bit and try again
					ThreadTools.sleep(100, TimeUnit.MILLISECONDS);
					continue;
				}

				process(nodeInfo);
			}
		}

		public void waitForFinish() {
			try {
				join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	public final TaskExecutor cpuTasks;
	public final SeqDB seqdb;
	public final NodeDB nodedb;
	public final StateInfo[] stateInfos;
	public final boolean includeStaticStatic;

	public final ConfEnergyCalculator[] ecalcs;

	private final List<MinimizationQueue> minimizationQueues;
	private final List<NodeThread> nodeThreads = new ArrayList<>();

	public NodeProcessor(TaskExecutor cpuTasks, SeqDB seqdb, NodeDB nodedb, StateInfo[] stateInfos, boolean includeStaticStatic, Parallelism parallelism, Structs.Precision precision) {

		this.cpuTasks = cpuTasks;
		this.seqdb = seqdb;
		this.nodedb = nodedb;
		this.stateInfos = stateInfos;
		this.includeStaticStatic = includeStaticStatic;

		// make the energy calculators
		ecalcs = Arrays.stream(stateInfos)
			.map(stateInfo -> ConfEnergyCalculator.makeBest(stateInfo.config.confSpace, parallelism, precision))
			.toArray(ConfEnergyCalculator[]::new);

		// set up the minimization queues
		minimizationQueues = Arrays.stream(stateInfos)
			.map(stateInfo -> new MinimizationQueue(stateInfo, ecalcs[stateInfo.config.state.index].maxBatchSize()))
			.collect(Collectors.toList());
	}

	@Override
	public void close() {

		// wait for the node threads to finish first
		for (var t : nodeThreads) {
			t.waitForFinish();
		}

		for (var ecalc : ecalcs) {
			ecalc.close();
		}
	}

	public void startNodeThreads(int numThreads, Directions directions) {
		for (int i=0; i<numThreads; i++) {
			nodeThreads.add(new NodeThread(i, directions));
		}
	}

	private void log(String msg, Object ... args) {
		nodedb.member.log(msg, args);
	}

	private Sequence makeSeq(int statei, int[] conf) {
		if (stateInfos[statei].config.state.isSequenced) {
			return seqdb.confSpace.seqSpace.makeSequence(stateInfos[statei].config.confSpace, conf);
		} else {
			return null;
		}
	}

	private enum Result {

		GotNode(true),
		NoState(false),
		NoTree(false),
		NoNode(false);

		public final boolean gotNode;

		Result(boolean gotNode) {
			this.gotNode = gotNode;
		}
	}

	private static class NodeInfo {

		final Result result;
		final NodeIndex.Node node;
		final RCs tree;

		NodeInfo(Result result, NodeIndex.Node node, RCs tree) {
			this.result = result;
			this.node = node;
			this.tree = tree;
		}

		NodeInfo(Result result) {
			this(result, null, null);
		}
	}

	private NodeInfo getNode(Directions directions) {

		// get the currently focused state
		int statei = directions.getFocusedStatei();
		if (statei < 0) {
			return new NodeInfo(Result.NoState);
		}

		// get the tree for this state
		RCs tree = directions.getTree(statei);
		if (tree == null) {
			return new NodeInfo(Result.NoTree);
		}

		// get the next node from that state
		var node = nodedb.removeHigh(statei);
		if (node == null) {
			return new NodeInfo(Result.NoNode);
		}

		// all is well
		return new NodeInfo(Result.GotNode, node, tree);
	}

	private NodeInfo process(NodeInfo nodeInfo) {

		if (nodeInfo.node.isLeaf()) {

			// see if this node's score is roughly as good as current predictions
			var currentScore = nodedb.perf.score(nodeInfo.node);
			if (nodeInfo.node.score.exp - currentScore.exp > 1) {

				// nope, it should have a much worse score
				// re-score it and put it back into nodedb
				nodedb.add(new NodeIndex.Node(nodeInfo.node, currentScore));

			} else {

				// the score looks good, minimize it
				minimize(nodeInfo.node);
			}
		} else {

			// interior node, expand it
			expand(nodeInfo);
		}

		return nodeInfo;
	}

	private void minimize(NodeIndex.Node node) {

		// the score look good, add it to the minimization queue
		var q = minimizationQueues.get(node.statei);
		List<NodeIndex.Node> nodes = q.add(node);

		// if we finished a batch, minimize it now
		if (nodes != null) {
			minimize(node.statei, nodes);
		}
	}

	private List<NodeIndex.Node> minimize(int statei, List<NodeIndex.Node> nodes) {

		// got a full batch, minimize it!
		var stateInfo = stateInfos[statei];

		// collect timing info for the minimizations
		Stopwatch stopwatch = new Stopwatch().start();

		// minimize the conformations
		var jobs = nodes.stream()
			.map(node -> new ConfEnergyCalculator.MinimizationJob(node.conf, makeInters(stateInfo, node.conf)))
			.collect(Collectors.toList());
		ecalcs[statei].minimizeEnergies(jobs);

		// TODO: if the GPUs are busy, minimize on the CPU

		// compute the bound energies for each conf
		double[] bounds = nodes.stream()
			.mapToDouble(node -> stateInfo.bcalc.freeEnergyPrecise(node.zSumUpper))
			.toArray();

		// update stats on energy bounds
		for (int i=0; i<nodes.size(); i++) {
			stateInfo.energyBoundStats.add(bounds[i], jobs.get(i).energy);
		}

		// update seqdb with boltzmann-weighted energies
		var seqdbBatch = seqdb.batch();
		for (int i=0; i<nodes.size(); i++) {
			var n = nodes.get(i);
			var e = jobs.get(i).energy;
			var z = stateInfo.bcalc.calcPrecise(e);
			seqdbBatch.addZConf(
				stateInfo.config.state,
				makeSeq(n.statei, n.conf),
				z,
				n.zSumUpper,
				new ConfSearch.EnergiedConf(n.conf, bounds[i], e)
			);
		}
		seqdbBatch.save();

		// update node performance
		stopwatch.stop();
		for (int i=0; i<nodes.size(); i++) {
			var n = nodes.get(i);
			var reduction = n.zSumUpper;
			long ns = stopwatch.getTimeNs()/nodes.size();
			nodedb.perf.update(n, ns, reduction);

			/* TEMP
			log("            minimization:   r %s  %10s  r/t %s",
				reduction,
				stopwatch.getTime(2),
				Log.formatBigEngineering(processor.seqdb.bigMath()
					.set(reduction)
					.div(stopwatch.getTimeNs())
					.get()
				)
			);
			*/
		}

		return nodes;
	}

	private List<PosInter> makeInters(StateInfo stateInfo, int[] conf) {
		if (includeStaticStatic) {
			return stateInfo.config.posInterGen.all(stateInfo.config.confSpace, conf);
		} else {
			return stateInfo.config.posInterGen.dynamic(stateInfo.config.confSpace, conf);
		}
	}

	public List<ConfEnergyCalculator.EnergiedCoords> minimizeCoords(int statei, List<int[]> confs) {

		var stateInfo = stateInfos[statei];

		var energiedCoords = confs.stream()
			.map(conf -> (ConfEnergyCalculator.EnergiedCoords)null)
			.collect(Collectors.toList());

		// minimize the confs
		// (can't use batches because we need the coords)
		for (int i=0; i<confs.size(); i++) {
			final int fi = i;
			cpuTasks.submit(
				() -> {
					int[] conf = confs.get(fi);
					var inters = makeInters(stateInfo, conf);
					energiedCoords.set(fi, ecalcs[statei].minimize(conf, inters));
					return 42;
				},
				answer -> {}
			);
		}

		cpuTasks.waitForFinish();

		return energiedCoords;
	}

	private void expand(NodeInfo nodeInfo) {

		// get timing info for the node expansion
		Stopwatch stopwatch = new Stopwatch().start();

		int statei = nodeInfo.node.statei;
		var stateInfo = stateInfos[statei];
		var confIndex = stateInfo.makeConfIndex();
		Conf.index(nodeInfo.node.conf, confIndex);

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
		int posi = stateInfo.posPermutation[confIndex.numDefined];

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

		seqdbBatch.save();

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
