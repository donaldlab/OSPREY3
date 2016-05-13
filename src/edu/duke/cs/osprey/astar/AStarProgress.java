package edu.duke.cs.osprey.astar;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;

import edu.duke.cs.osprey.tools.Stopwatch;

public class AStarProgress {
	
	private static final int ReportIntervalMs = 1 * 1000; // TEMP
	
	private int numLevels;
	private int level;
	private int deepestLevel;
	private double score;
	private long numNodesInQueue;
	private long numNodesExplored;
	private long numNodesExploredThisReport;
	private Stopwatch stopwatch;
	private int msRunning;
	
	public void init() {
		numLevels = 0;
		level = 0;
		deepestLevel = 0;
		score = Double.POSITIVE_INFINITY;
		numNodesInQueue = 0;
		numNodesExplored = 0;
		numNodesExploredThisReport = 0;
		stopwatch = new Stopwatch();
		msRunning = 0;
	}
	
	public void reportNode(AStarNode node, int numNodesInQueue) {
		
		this.numLevels = node.getNodeAssignments().length;
		this.level = getLevel(node.getNodeAssignments());
		this.deepestLevel = Math.max(this.deepestLevel, this.level);
		this.score = node.getScore();
		this.numNodesInQueue = numNodesInQueue;
		this.numNodesExplored++;
		this.numNodesExploredThisReport++;
	
		if (!stopwatch.isRunning()) {
			stopwatch.start();
			msRunning = 0;
		}
		
		// should we write a progress report?
		int msRunning = (int)stopwatch.getTimeMs();
		if (msRunning >= this.msRunning + ReportIntervalMs) {
			printProgressReport();
			this.msRunning = msRunning;
			this.numNodesExploredThisReport = 0;
		}
	}
	
	public void printProgressReport() {
		// TODO: configurable output? logging framework?
		System.out.println(makeProgressReport());
	}
	
	public String makeProgressReport() {
		double diffMs = stopwatch.getTimeMs() - this.msRunning;
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		return String.format("AStar Progress: score=%18.12f, level=%4d/%4d/%4d, queuedNodes=%14d, exploredNodes=%14d, nodes/sec=%5d, time=%s, heapMem=%.0f%%",
			score,
			level, deepestLevel, numLevels,
			numNodesInQueue, numNodesExplored,
			(int)(numNodesExploredThisReport*1000/diffMs),
			stopwatch.getTime(2),
			100f*heapMem.getUsed()/heapMem.getMax()
		);
	}
	
	private int getLevel(int[] conf) {
		int level = 0;
		for (int i=0; i<conf.length; i++) {
			if (conf[i] >= 0) {
				
				// defined conf at this position, increment level
				level++;
			}
		}
		return level;
	}
}
