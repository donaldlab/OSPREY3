package edu.duke.cs.osprey.astar;

import java.io.Serializable;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;

import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.tools.JvmMem;
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
	private int numLeafNodes;
	private double goalScore;
	
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
		numLeafNodes = 0;
	}
	
	public void reportInternalNode(int level, double gscore, double hscore, long numNodesInQueue, int numAddedToQueue) {
		
		this.numNodesExpanded++;
		this.numNodesInQueue = numNodesInQueue;
		this.numNodesQueuedThisReport += numAddedToQueue;
	
		// if we've hit a leaf node, only track node expansion stats
		if (numLeafNodes > 0) {
			return;
		}
		
		this.level = level;
		this.deepestLevel = Math.max(this.deepestLevel, this.level);
		this.gscore = gscore;
		this.hscore = hscore;
		
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
		return String.format("A* g:%10.4f, h:%10.4f, f:%10.4f, level:%4d/%4d/%4d, expanded:%10d, queued:%10d, scored/sec:%5d, time:%s, heapMem:%s, extMem:%s",
			gscore, hscore, gscore + hscore,
			level, deepestLevel, numLevels - 1,
			numNodesExpanded, numNodesInQueue,
			(int)(numNodesQueuedThisReport*1000/diffMs),
			stopwatch.getTime(2),
			JvmMem.getOldPool(),
			ExternalMemory.getUsageReport()
		);
	}
	
	public void setGoalScore(double val) {
		goalScore = val;
	}
	
	public void reportLeafNode(double gscore, long numNodesInQueue) {
		
		this.numNodesInQueue = numNodesInQueue;
		
		// if this is the first leaf node, print one last progress report
		if (this.numLeafNodes == 0) {
			
			this.level = numLevels - 1;
			this.deepestLevel = numLevels - 1;
			this.gscore = gscore;
			this.hscore = 0;
			
			printProgressReport();
		}
		
		this.gscore = gscore;
		this.numLeafNodes++;
		
		// should we write a progress report?
		int msRunning = (int)stopwatch.getTimeMs();
		if (msRunning >= this.msRunning + ReportIntervalMs) {
			printLeafProgressReport();
			this.msRunning = msRunning;
			this.numNodesQueuedThisReport = 0;
		}
	}
	
	private void printLeafProgressReport() {
		// TODO: configurable output? logging framework?
		System.out.println(makeLeafProgressReport());
	}

	private String makeLeafProgressReport() {
		double diffMs = stopwatch.getTimeMs() - this.msRunning;
		return String.format("A* leaf nodes:%10d, score:%14.8f, remaining:%14.8f, expanded:%10d, queued:%10d, scored/sec:%5d, time:%s, heapMem:%s, extMem:%s",
			numLeafNodes,
			gscore, goalScore - gscore,
			numNodesExpanded, numNodesInQueue,
			(int)(numNodesQueuedThisReport*1000/diffMs),
			stopwatch.getTime(2),
			JvmMem.getOldPool(),
			ExternalMemory.getUsageReport()
		);
	}
}
