package edu.duke.cs.osprey.astar;

import java.io.Serializable;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;

import edu.duke.cs.osprey.tools.Stopwatch;

public class AStarProgress implements Serializable {

	private static final long serialVersionUID = 767988861537797830L;

	private static final int ReportIntervalMs = 10 * 1000; // TODO: make configurable
	
	private int numLevels;
	private int level;
	private int deepestLevel;
	private double gscore;
	private double hscore;
	private long numNodesExpanded;
	private long numNodesInQueue;
	private long numNodesQueuedThisReport;
	private Stopwatch stopwatch;
	private int msRunning;
	
	public AStarProgress(int numLevels) {
		this.numLevels = numLevels;
		level = 0;
		deepestLevel = 0;
		gscore = Double.POSITIVE_INFINITY;
		hscore = Double.POSITIVE_INFINITY;
		numNodesInQueue = 0;
		numNodesQueuedThisReport = 0;
		stopwatch = new Stopwatch();
		msRunning = 0;
	}
	
	public void reportNode(int level, double gscore, double hscore, int numNodesInQueue, int numAddedToQueue) {
		
		this.level = level;
		this.deepestLevel = Math.max(this.deepestLevel, this.level);
		this.gscore = gscore;
		this.hscore = hscore;
		this.numNodesExpanded++;
		this.numNodesInQueue = numNodesInQueue;
		this.numNodesQueuedThisReport += numAddedToQueue;
	
		if (!stopwatch.isRunning()) {
			stopwatch.start();
			msRunning = 0;
		}
		
		// should we write a progress report?
		int msRunning = (int)stopwatch.getTimeMs();
		if (msRunning >= this.msRunning + ReportIntervalMs) {
			printProgressReport();
			this.msRunning = msRunning;
			this.numNodesQueuedThisReport = 0;
		}
	}
	
	public void printProgressReport() {
		// TODO: configurable output? logging framework?
		System.out.println(makeProgressReport());
	}
	
	public String makeProgressReport() {
		double diffMs = stopwatch.getTimeMs() - this.msRunning;
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		return String.format("A* gscore:%14.8f, hscore:%14.8f, level:%4d/%4d/%4d, expanded:%10d, queued:%10d, scored/sec:%5d, time:%s, heapMem:%.0f%%",
			gscore, hscore,
			level, deepestLevel, numLevels - 1,
			numNodesExpanded, numNodesInQueue,
			(int)(numNodesQueuedThisReport*1000/diffMs),
			stopwatch.getTime(2),
			100f*heapMem.getUsed()/heapMem.getMax()
		);
	}
}
