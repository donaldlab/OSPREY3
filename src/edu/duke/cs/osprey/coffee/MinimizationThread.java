package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class MinimizationThread {

	public final StateInfo stateInfo;
	public final NodeProcessor processor;
	public final TaskExecutor tasks;

	private final int batchSize;
	private final int queueSize;
	private final Deque<NodeIndex.Node> queue;
	private final AtomicBoolean isRunning;

	private final Thread thread;

	public MinimizationThread(StateInfo stateInfo, NodeProcessor processor, TaskExecutor tasks) {

		this.stateInfo = stateInfo;
		this.processor = processor;
		this.tasks = tasks;

		batchSize = stateInfo.config.ecalc.maxBatchSize();
		queueSize = batchSize*2;
		queue = new ArrayDeque<>(queueSize);
		isRunning = new AtomicBoolean(true);

		thread = new Thread(() -> run());
		thread.setName("MinimizationThread-" + stateInfo.config.state.name);
		thread.setDaemon(false);
		thread.start();
	}

	public void askToStop() {
		isRunning.set(false);
	}

	public void join() {
		try {
			thread.join();
		} catch (InterruptedException ex) {
			throw new RuntimeException("can't join thread", ex);
		}
	}

	public void addNode(NodeIndex.Node node) {
		while (true) {

			synchronized (queue) {
				if (queue.size() < queueSize) {
					queue.addLast(node);
					break;
				}
			}

			// queue is full, wait a bit and try again
			ThreadTools.sleep(100, TimeUnit.MILLISECONDS);
		}
	}

	private void run() {

		while (isRunning.get()) {

			// get a batch of minimizations
			List<NodeIndex.Node> nodeBatch;
			synchronized (queue) {
				if (queue.isEmpty()) {
					nodeBatch = null;
				} else {
					nodeBatch = new ArrayList<>(batchSize);
					for (int i=0; i<batchSize; i++) {
						NodeIndex.Node node = queue.pollFirst();
						if (node == null) {
							break;
						}
						nodeBatch.add(node);
					}
				}
			}

			// no nodes? wait a bit and try again
			if (nodeBatch == null) {
				ThreadTools.sleep(100, TimeUnit.MILLISECONDS);
				continue;
			}

			// minimize the nodes
			tasks.submit(
				() -> {

					// collect timing info for the minimizations
					Stopwatch stopwatch = new Stopwatch().start();

					// actually do the minimizations
					var zs = stateInfo.zPaths(nodeBatch.stream()
						.map(node -> node.conf)
						.collect(Collectors.toList()));

					// prep the seqdb batch
					var seqdbBatch = processor.seqdb.batch();
					for (int i=0; i<nodeBatch.size(); i++) {
						var node = nodeBatch.get(i);
						seqdbBatch.addZConf(
							stateInfo.config.state,
							processor.makeSeq(node.statei, node.conf),
							zs.get(i),
							node.zSumUpper
						);
					}

					// update node performance
					stopwatch.stop();
					for (int i=0; i<nodeBatch.size(); i++) {
						var node = nodeBatch.get(i);
						var reduction = new BigExp(node.zSumUpper);
						reduction.sub(zs.get(i));
						long ns = stopwatch.getTimeNs()/nodeBatch.size();
						processor.nodedb.perf.update(node, ns, reduction);

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

					seqdbBatch.save();

					return 42;
				},
				theAnswer -> {}
			);
		}
	}
}
