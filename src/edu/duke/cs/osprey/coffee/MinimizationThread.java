package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.coffee.nodedb.NodeIndex;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadTools;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.Collectors;


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

			// minimize them
			tasks.submit(
				() -> stateInfo.zPaths(nodeBatch.stream()
					.map(node -> node.conf)
					.collect(Collectors.toList())),
				zs -> {

					// update seqdb
					var seqdbBatch = processor.seqdb.batch();
					for (int i=0; i<nodeBatch.size(); i++) {
						var node = nodeBatch.get(i);
						seqdbBatch.addZConf(
							stateInfo.config.state,
							processor.makeSeq(node.statei, node.conf),
							zs.get(i),
							node.score
						);
					}
					seqdbBatch.save();
				}
			);
		}
	}
}
