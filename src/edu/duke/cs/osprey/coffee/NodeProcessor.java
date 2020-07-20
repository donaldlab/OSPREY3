package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeDB;
import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.coffee.seqdb.Batch;
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

	private static class FlushTracker {

		final long flushNs = TimeUnit.MILLISECONDS.toNanos(100);

		int lastStatei = -1;
		long lastFlush = -1;

		boolean stateChanged(int statei) {
			return statei != lastStatei;
		}

		void setState(int statei) {
			lastStatei = statei;
		}

		boolean shouldFlush() {
			return System.nanoTime() >= lastFlush + flushNs;
		}

		void flushed() {
			lastFlush = System.nanoTime();
		}
	}

	private class NodeThread extends Thread {

		final int id;
		final Directions directions;

		Batch seqBatch = null;
		final List<NodeIndex.Node> nodesIncoming = new ArrayList<>();
		final List<NodeIndex.Node> nodesOutgoing = new ArrayList<>();
		final FlushTracker flushTracker = new FlushTracker();

		final int nodeBatchSize = 100;

		NodeThread(int id, Directions directions) {

			this.id = id;
			this.directions = directions;

			setName("Node-" + id);
			setDaemon(true);
			start();
		}

		@Override
		public void run() {

			Runnable waitABit = () -> ThreadTools.sleep(100, TimeUnit.MILLISECONDS);

			if (seqdb != null) {
				seqBatch = seqdb.batch();
			}

			while (directions.isRunning()) {

				// flush if needed
				if (flushTracker.shouldFlush()) {
					flush();
				}

				long startNs = System.nanoTime();

				// get the currently focused state
				int statei = directions.getFocusedStatei();
				if (statei < 0) {
					waitABit.run();
					continue;
				}

				// if the state changed, flush
				if (flushTracker.stateChanged(statei)) {
					flush();
					flushTracker.setState(statei);
				}

				// get the tree for this state
				RCs tree = directions.getTree(statei);
				if (tree == null) {
					waitABit.run();
					continue;
				}

				// get the next nodes from that state
				nodedb.removeHigh(statei, nodeBatchSize, nodesIncoming);
				if (nodesIncoming.isEmpty()) {
					waitABit.run();
					continue;
				}

				// got some nodes! process them!
				long nodeNs = System.nanoTime() - startNs;
				for (var node : nodesIncoming) {

					process(new NodeInfo(node, tree, nodeNs/nodesIncoming.size()), seqBatch, nodesOutgoing);

					if (flushTracker.shouldFlush()) {
						flush();
					}
				}
				nodesIncoming.clear();
			}
		}

		private void flush() {

			if (seqBatch != null) {
				seqBatch.save();
			}

			if (!nodesOutgoing.isEmpty()) {
				assert (flushTracker.lastStatei >= 0);
				nodedb.add(flushTracker.lastStatei, nodesOutgoing);
				nodesOutgoing.clear();
			}

			flushTracker.flushed();
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

		Batch seqBatch = null;
		FlushTracker flushTracker = new FlushTracker();

		GpuThread(int id, Directions directions) {

			this.id = id;
			this.directions = directions;

			setName("GpuMinimizer-" + id);
			setDaemon(true);
			start();
		}

		@Override
		public void run() {

			if (seqdb != null) {
				seqBatch = seqdb.batch();
			}

			while (directions.isRunning()) {

				// flush if we haven't done it in a while
				if (flushTracker.shouldFlush()) {
					flush();
				}

				// listen to directions
				int statei = directions.getFocusedStatei();
				if (statei < 0) {
					ThreadTools.sleep(500, TimeUnit.MILLISECONDS);
					continue;
				}

				// if the state changed, flush
				if (flushTracker.stateChanged(statei)) {
					flush();
					flushTracker.setState(statei);
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

		void waitForFinish() {
			try {
				join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}

		void minimize(StateInfo stateInfo, CudaConfEnergyCalculator ecalc, List<NodeInfo> nodes) {

			// collect timing info for the minimizations
			Stopwatch stopwatch = new Stopwatch().start();

			// minimize the nodes
			var jobs = nodes.stream()
				.map(info -> new ConfEnergyCalculator.MinimizationJob(info.node.conf, makeInters(stateInfo, info.node.conf)))
				.collect(Collectors.toList());
			ecalc.minimizeEnergies(jobs);

			minimized(stateInfo, nodes, jobs, stopwatch, seqBatch);
		}

		void flush() {

			if (seqBatch != null) {
				seqBatch.save();
			}

			flushTracker.flushed();
		}
	}

	private static class MinimizationQueue {

		public final int capacity;
		public final int batchSize;

		private final Deque<NodeInfo> nodes;
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

		NodeInfo offer(NodeInfo node) {
			final ReentrantLock lock = this.lock;
			lock.lock();
			try {

				// if the queue is full, return the first node
				NodeInfo out = null;
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

		List<NodeInfo> poll(int count, long timeout, TimeUnit unit) {
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
					var batch = new ArrayList<NodeInfo>(count);
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

	private class DropThread extends Thread {

		final Directions directions;

		final Deque<NodeIndex.Node> nodeQueue;
		final ReentrantLock lock;
		final Condition nodesReady;

		Batch seqBatch = null;
		FlushTracker flushTracker = new FlushTracker();

		public DropThread(Directions directions) {

			this.directions = directions;

			nodeQueue = new ArrayDeque<>();
			lock = new ReentrantLock(false);
			nodesReady = lock.newCondition();

			setName("NodeDrops");
			setDaemon(true);
			start();
		}

		@Override
		public void run() {

			if (seqdb != null) {
				seqBatch = seqdb.batch();
			}

			while (directions.isRunning()) {

				// flush if we haven't done it in a while
				if (flushTracker.shouldFlush()) {
					flush();
				}

				// listen to directions
				int statei = directions.getFocusedStatei();
				if (statei < 0) {
					ThreadTools.sleep(500, TimeUnit.MILLISECONDS);
					continue;
				}

				// if the state changed, flush
				if (flushTracker.stateChanged(statei)) {
					flush();
					flushTracker.setState(statei);
				}

				// check the drop queue
				try {
					long timeoutNs = TimeUnit.MILLISECONDS.toNanos(100);
					final ReentrantLock lock = this.lock;
					lock.lockInterruptibly();
					try {

						// wait for nodes to drop, if any
						while (nodeQueue.isEmpty()) {
							if (timeoutNs <= 0L) {
								break;
							}
							timeoutNs = nodesReady.awaitNanos(timeoutNs);
						}

						// drop any nodes we got
						while (!nodeQueue.isEmpty()) {
							var node = nodeQueue.poll();
							seqBatch.drop(
								stateInfos[node.statei].config.state,
								makeSeq(node.statei, node.conf),
								node.zSumUpper
							);
						}

					} finally {
						lock.unlock();
					}
				} catch (InterruptedException ex) {
					throw new RuntimeException(ex);
				}
			}
		}

		void waitForFinish() {
			try {
				join();
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}

		void addNodes(Stream<NodeIndex.Node> nodes) {
			final ReentrantLock lock = this.lock;
			lock.lock();
			try {
				nodes.forEach(node -> nodeQueue.offer(node));

				// signal pollers if needed
				if (!nodeQueue.isEmpty()) {
					nodesReady.signal();
				}

			} finally {
				lock.unlock();
			}
		}

		void flush() {
			if (seqBatch != null) {
				seqBatch.save();
			}
			flushTracker.flushed();
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
	private DropThread dropThread = null;

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
		dropThread.waitForFinish();

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

		// start the drop thread
		dropThread = new DropThread(directions);
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

		final NodeIndex.Node node;
		final RCs tree;
		final long aquisitionNs;

		NodeInfo(NodeIndex.Node node, RCs tree, long aquisitionNs) {
			this.node = node;
			this.tree = tree;
			this.aquisitionNs = aquisitionNs;
		}
	}

	private void process(NodeInfo nodeInfo, Batch seqBatch, List<NodeIndex.Node> nodeBatch) {

		if (nodeInfo.node.isLeaf()) {

			// see if this node's score is roughly as good as current predictions
			var currentScore = nodedb.perf.score(nodeInfo.node);
			if (nodeInfo.node.score.exp - currentScore.exp > 1) {

				// nope, it should have a much worse score
				// re-score it and put it back into nodedb
				nodeBatch.add(new NodeIndex.Node(nodeInfo.node, currentScore));

			} else {

				// the score looks good, minimize it
				minimize(nodeInfo, seqBatch);
			}
		} else {

			// interior node, expand it
			expand(nodeInfo, seqBatch, nodeBatch);
		}
	}

	private void minimize(NodeInfo nodeInfo, Batch seqBatch) {

		int statei = nodeInfo.node.statei;
		var stateInfo = stateInfos[nodeInfo.node.statei];

		// do we have GPUs?
		if (gpuEcalcs != null) {

			// yup, put the node on the queue and let the GPUs deal with it
			nodeInfo = minimizationQueues.get(statei).offer(nodeInfo);
			if (nodeInfo == null) {
				return;
			}
		}

		// we don't have GPUs or they're busy, so minimize on this CPU thread

		// collect timing info for the minimizations
		Stopwatch stopwatch = new Stopwatch().start();

		// minimize it
		var nodeInfos = Collections.singletonList(nodeInfo);
		var jobs = nodeInfos.stream()
			.map(info -> new ConfEnergyCalculator.MinimizationJob(info.node.conf, makeInters(stateInfo, info.node.conf)))
			.collect(Collectors.toList());

		cpuEcalcs[statei].minimizeEnergies(jobs);

		minimized(stateInfo, nodeInfos, jobs, stopwatch, seqBatch);
	}

	private List<PosInter> makeInters(StateInfo stateInfo, int[] conf) {
		return stateInfo.config.makeInters(conf, includeStaticStatic);
	}

	private void minimized(StateInfo stateInfo, List<NodeInfo> nodeInfos, List<ConfEnergyCalculator.MinimizationJob> jobs, Stopwatch stopwatch, Batch seqBatch) {

		// compute the bound energies for each conf
		double[] bounds = nodeInfos.stream()
			.mapToDouble(info -> stateInfo.zmat.bcalc.freeEnergyPrecise(info.node.zSumUpper))
			.toArray();

		// update stats on energy bounds
		for (int i=0; i<nodeInfos.size(); i++) {
			stateInfo.energyBoundStats.add(bounds[i], jobs.get(i).energy);
		}

		// update seqdb with boltzmann-weighted energies if needed
		if (seqdb != null) {
			for (int i=0; i<nodeInfos.size(); i++) {
				var n = nodeInfos.get(i).node;
				var e = jobs.get(i).energy;
				var z = stateInfo.zmat.bcalc.calcPrecise(e);
				seqBatch.addZConf(
					stateInfo.config.state,
					makeSeq(n.statei, n.conf),
					z,
					n.zSumUpper,
					new ConfSearch.EnergiedConf(n.conf, bounds[i], e)
				);
			}
		}

		// update node performance
		stopwatch.stop();
		for (var nodeInfo : nodeInfos) {
			var reduction = nodeInfo.node.zSumUpper;
			long ns = stopwatch.getTimeNs()/nodeInfos.size() + nodeInfo.aquisitionNs;
			nodedb.perf.update(nodeInfo.node, ns, reduction);
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

	public void expand(NodeIndex.Node node, RCs tree, List<NodeIndex.Node> nodeBatch) {
		var nodeInfo = new NodeInfo(node, tree, 0);
		expand(nodeInfo, null, nodeBatch);
	}

	private void expand(NodeInfo nodeInfo, Batch seqBatch, List<NodeIndex.Node> nodeBatch) {

		// just in case ...
		if (nodeInfo.node.isLeaf()) {
			throw new IllegalArgumentException("can't expand leaf node");
		}

		// get timing info for the node expansion
		Stopwatch stopwatch = new Stopwatch().start();

		int statei = nodeInfo.node.statei;
		var stateInfo = stateInfos[statei];
		var confIndex = stateInfo.makeConfIndex();
		Conf.index(nodeInfo.node.conf, confIndex);

		// remove the old bound from SeqDB if needed
		if (seqBatch != null) {
			seqBatch.subZSumUpper(
				stateInfo.config.state,
				makeSeq(statei, nodeInfo.node.conf),
				nodeInfo.node.zSumUpper
			);
		}

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
			nodeBatch.add(new NodeIndex.Node(
				statei, conf, zSumUpper,
				nodedb.perf.score(statei, conf, zSumUpper)
			));

			// add the new bound if needed
			if (seqBatch != null) {
				seqBatch.addZSumUpper(
					stateInfo.config.state,
					makeSeq(statei, conf),
					zSumUpper
				);
			}

			// update the reduction calculation
			reduction.sub(zSumUpper);

			confIndex.unassignInPlace(posi);
		}

		// update the node performance
		stopwatch.stop();
		nodedb.perf.update(nodeInfo.node, stopwatch.getTimeNs() + nodeInfo.aquisitionNs, reduction);
	}

	public void handleDrops(Stream<NodeIndex.Node> nodes) {

		// NOTE: this is called on the NodeDB thread!

		if (seqdb != null) {
			// move the nodes to the drop thread
			dropThread.addNodes(nodes);
		}
	}
}
