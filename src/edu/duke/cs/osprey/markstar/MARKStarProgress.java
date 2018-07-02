package edu.duke.cs.osprey.markstar;

import java.io.Serializable;

import edu.duke.cs.osprey.astar.AStarProgress;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.tools.JvmMem;
import edu.duke.cs.osprey.tools.Stopwatch;

public class MARKStarProgress extends AStarProgress {

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
    private double curEpsilon;
    private int numPartialMinimizations;

    public MARKStarProgress(int numLevels) {
        super(numLevels);
        this.numLevels = numLevels;
        level = 0;
        deepestLevel = 0;
        gscore = Double.POSITIVE_INFINITY;
        hscore = Double.POSITIVE_INFINITY;
        numNodesInQueue = 0;
        numNodesQueuedThisReport = 0;
        numPartialMinimizations = 0;
        stopwatch = new Stopwatch();
        msRunning = 0;
        numLeafNodes = 0;
        curEpsilon = 1;
    }

    public void reportInternalNode(int level, double gscore, double hscore, long numNodesInQueue, int numAddedToQueue, double curEpsilon) {

        this.numNodesExpanded++;
        this.numNodesInQueue = numNodesInQueue;
        this.numNodesQueuedThisReport += numAddedToQueue;
        this.curEpsilon = curEpsilon;

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

    public String makeProgressReport() {
        double diffMs = stopwatch.getTimeMs() - this.msRunning;
        return String.format("MARK* expanded:%10d,"
                +" queued:%10d, scored/sec:%5d, partial minimizations:%5d, time:%s, delta:%10.4f heapMem:%s, extMem:%s",
                numNodesExpanded, numNodesInQueue,
                (int)(numNodesQueuedThisReport*1000/diffMs),
                numPartialMinimizations,
                stopwatch.getTime(2),
                curEpsilon,
                JvmMem.getOldPool(),
                ExternalMemory.getUsageReport()
        );
    }

    public void setGoalScore(double val) {
        goalScore = val;
    }

    public void reportLeafNode(double gscore, long numNodesInQueue, double curEpsilon) {
        this.curEpsilon = curEpsilon;
        super.reportLeafNode(gscore, numNodesInQueue);
    }

    protected String makeLeafProgressReport() {
        double diffMs = stopwatch.getTimeMs() - this.msRunning;
        return String.format("A* leaf nodes:%10d, score:%14.8f, remaining:%14.8f, expanded:%10d, queued:%10d,"
                +" scored/sec:%5d, time:%s, delta:%10.4f, heapMem:%s, extMem:%s",
                numLeafNodes,
                gscore, goalScore - gscore,
                numNodesExpanded, numNodesInQueue,
                (int)(numNodesQueuedThisReport*1000/diffMs),
                stopwatch.getTime(2),
                curEpsilon,
                JvmMem.getOldPool(),
                ExternalMemory.getUsageReport()
        );
    }

    public void reportPartialMinimization(int numPartialMinimizations, double curEpsilon) {
        this.numPartialMinimizations += numPartialMinimizations;

        this.curEpsilon = curEpsilon;

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
            this.numPartialMinimizations = 0;
        }

    }
}
