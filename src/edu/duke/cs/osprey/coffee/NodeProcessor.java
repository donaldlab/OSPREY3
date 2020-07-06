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
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.stream.Stream;


public class NodeProcessor implements AutoCloseable {

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

	private class GpuThread extends Thread {

		final int id;
		final Directions directions;

		GpuThread(int id, Directions directions) {

			this.id = id;
			this.directions = directions;

			setName("GpuMinimizer-" + id);
			setDaemon(true);
			start();
		}

		@Override
		public void run() {

			while (directions.isRunning()) {

				// listen to directions
				int statei = directions.getFocusedStatei();
				if (statei < 0) {
					ThreadTools.sleep(500, TimeUnit.MILLISECONDS);
					continue;
				}

				var stateInfo = stateInfos[statei];
				var q = minimizationQueues.get(statei);
				var ecalc = gpuEcalcs[statei];

				// get the next batch to minimize
				var nodes = q.poll(ecalc.maxBatchSize(), 100, TimeUnit.MILLISECONDS);
				if (nodes != null) {
					minimize(stateInfo, ecalc, nodes);
				}
			}
		}

		public void waitForFinish() {
			try {
				join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}

		private void minimize(StateInfo stateInfo, CudaConfEnergyCalculator ecalc, List<NodeIndex.Node> nodes) {

			// collect timing info for the minimizations
			Stopwatch stopwatch = new Stopwatch().start();

			// minimize the nodes
			var jobs = nodes.stream()
				.map(n -> new ConfEnergyCalculator.MinimizationJob(n.conf, makeInters(stateInfo, n.conf)))
				.collect(Collectors.toList());
			ecalc.minimizeEnergies(jobs);

			minimized(stateInfo, nodes, jobs, stopwatch);
		}
	}

	private static class MinimizationQueue {

		public final int capacity;
		public final int batchSize;

		private final Deque<NodeIndex.Node> nodes;
		private final ReentrantLock lock;
		private final Condition batchReady;

		MinimizationQueue(int capacity, int batchSize) {

			this.capacity = capacity;
			this.batchSize = batchSize;

			nodes = new ArrayDeque<>(capacity);
			lock = new ReentrantLock(false);
			batchReady = lock.newCondition();
		}

		int size() {
			final ReentrantLock lock = this.lock;
			lock.lock();
			try {
				return nodes.size();
			} finally {
				lock.unlock();
			}
		}

		NodeIndex.Node offer(NodeIndex.Node node) {
			final ReentrantLock lock = this.lock;
			lock.lock();
			try {

				// if the queue is full, return the first node
				NodeIndex.Node out = null;
				if (nodes.size() == capacity) {
					out = nodes.poll();
				}

				// add the node to the queue
				nodes.offer(node);

				// signal pollers if needed
				if (nodes.size() >= batchSize) {
					batchReady.signal();
				}

				return out;

			} finally {
				lock.unlock();
			}
		}

		List<NodeIndex.Node> poll(int count, long timeout, TimeUnit unit) {
			try {
				long timeoutNs = unit.toNanos(timeout);
				final ReentrantLock lock = this.lock;
				lock.lockInterruptibly();
				try {

					// wait for the batch to fill up
					while (nodes.size() < count) {
						if (timeoutNs <= 0L) {

							// no batch was ready in time
							return null;
						}
						timeoutNs = batchReady.awaitNanos(timeoutNs);
					}

					// poll off the batch
					var batch = new ArrayList<NodeIndex.Node>(count);
					while (batch.size() < count) {
						batch.add(nodes.poll());
					}
					return batch;

				} finally {
					lock.unlock();
				}
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
	public final Parallelism parallelism;

	public final ConfEnergyCalculator[] cpuEcalcs;
	public final CudaConfEnergyCalculator[] gpuEcalcs;

	private final List<NodeThread> nodeThreads = new ArrayList<>();
	private final List<GpuThread> gpuThreads = new ArrayList<>();
	private final List<MinimizationQueue> minimizationQueues = new ArrayList<>();

	public NodeProcessor(TaskExecutor cpuTasks, SeqDB seqdb, NodeDB nodedb, StateInfo[] stateInfos, boolean includeStaticStatic, Parallelism parallelism, Structs.Precision precision) {

		this.cpuTasks = cpuTasks;
		this.seqdb = seqdb;
		this.nodedb = nodedb;
		this.stateInfos = stateInfos;
		this.includeStaticStatic = includeStaticStatic;
		this.parallelism = parallelism;

		// make the energy calculators
		cpuEcalcs = Arrays.stream(stateInfos)
			.map(stateInfo -> new NativeConfEnergyCalculator(stateInfo.config.confSpace, precision))
			.toArray(ConfEnergyCalculator[]::new);
		if (parallelism.numGpus > 0) {
			gpuEcalcs = Arrays.stream(stateInfos)
				.map(stateInfo -> new CudaConfEnergyCalculator(stateInfo.config.confSpace, precision, parallelism))
				.toArray(CudaConfEnergyCalculator[]::new);
		} else {
			gpuEcalcs = null;
		}
	}

	@Override
	public void close() {

		// wait for the node threads to finish first
		for (var t : nodeThreads) {
			t.waitForFinish();
		}
		for (var t : gpuThreads) {
			t.waitForFinish();
		}

		// cleanup the ecalcs
		for (var ecalc : cpuEcalcs) {
			ecalc.close();
		}
		if (gpuEcalcs != null) {
			for (var ecalc : gpuEcalcs) {
				ecalc.close();
			}
		}
	}

	public void start(int numThreads, Directions directions) {

		if (!nodeThreads.isEmpty() || !gpuThreads.isEmpty()) {
			throw new IllegalStateException("threads already started");
		}

		// start the node threads
		for (int i=0; i<numThreads; i++) {
			nodeThreads.add(new NodeThread(i, directions));
		}

		// start the GPU threads too, if needed
		if (gpuEcalcs != null) {

			// all states should have the same GPU settings
			int numStreams = gpuEcalcs[0].numStreams();
			int batchSize = gpuEcalcs[0].maxBatchSize();

			// make the queues
			// make them big enough so all the GPU threads can get more work without waiting
			int queueCapacity = numStreams*batchSize*6;
			for (var ignored : stateInfos) {
				minimizationQueues.add(new MinimizationQueue(queueCapacity, batchSize));
			}

			// start the threads
			for (int streami=0; streami<numStreams; streami++) {
				gpuThreads.add(new GpuThread(streami, directions));
			}
		}
	}

	public int getMinimizationQueueSize(int statei) {
		if (minimizationQueues.isEmpty()) {
			return -1;
		}
		return minimizationQueues.get(statei).size();
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

		var stateInfo = stateInfos[node.statei];

		// do we have GPUs?
		if (gpuEcalcs != null) {

			// yup, put the node on the queue and let the GPUs deal with it
			node = minimizationQueues.get(node.statei).offer(node);
			if (node == null) {
				return;
			}
		}

		// we don't have GPUs or they're busy, so minimize on this CPU thread

		// collect timing info for the minimizations
		Stopwatch stopwatch = new Stopwatch().start();

		// minimize it
		var nodes = Collections.singletonList(node);
		var jobs = nodes.stream()
			.map(n -> new ConfEnergyCalculator.MinimizationJob(n.conf, makeInters(stateInfo, n.conf)))
			.collect(Collectors.toList());

		cpuEcalcs[node.statei].minimizeEnergies(jobs);

		minimized(stateInfo, nodes, jobs, stopwatch);
	}

	private List<PosInter> makeInters(StateInfo stateInfo, int[] conf) {
		if (includeStaticStatic) {
			return stateInfo.config.posInterGen.all(stateInfo.config.confSpace, conf);
		} else {
			return stateInfo.config.posInterGen.dynamic(stateInfo.config.confSpace, conf);
		}
	}

	private void minimized(StateInfo stateInfo, List<NodeIndex.Node> nodes, List<ConfEnergyCalculator.MinimizationJob> jobs, Stopwatch stopwatch) {

		// compute the bound energies for each conf
		double[] bounds = nodes.stream()
			.mapToDouble(n -> stateInfo.bcalc.freeEnergyPrecise(n.zSumUpper))
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
			var z = stateInfo.bcalc.calcPrecise(e); // TODO: is this bottlenecking the GPU threads?
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
	}

	public List<ConfEnergyCalculator.EnergiedCoords> minimizeCoords(int statei, List<int[]> confs) {

		var stateInfo = stateInfos[statei];

		var energiedCoords = confs.stream()
			.map(conf -> (ConfEnergyCalculator.EnergiedCoords)null)
			.collect(Collectors.toList());

		// pick an ecalc
		ConfEnergyCalculator ecalc;
		if (gpuEcalcs != null) {
			ecalc = gpuEcalcs[statei];
		} else {
			ecalc = cpuEcalcs[statei];
		}

		// minimize the confs
		// (can't use batches because we need the coords)
		for (int i=0; i<confs.size(); i++) {
			final int fi = i;
			cpuTasks.submit(
				() -> {
					int[] conf = confs.get(fi);
					var inters = makeInters(stateInfo, conf);
					energiedCoords.set(fi, ecalc.minimize(conf, inters));
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
