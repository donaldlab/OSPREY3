package edu.duke.cs.osprey.tools;

import edu.duke.cs.osprey.externalMemory.ExternalMemory;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;

public class Progress {
	
	private static final long NSpMS = 1000000;
	private static final int UpdateIntervalMs = 500; // half a second
	private static final int MaxLogSize = 120; // 60 seconds of log history
	
	private static final Model DefaultModel = Model.Linear;
	private static final int DefaultReportIntervalMs = 5*1000; // 5 seconds
	
	public static enum Model {
		
		Linear {
			@Override
			public long getEtaMs(Deque<Long> workLog, Deque<Long> timeMsLog, long totalWork) {
				
				// not enough data?
				if (workLog.size() < 2) {
					return -1;
				}
				
				// perform simple linear regression
				double sumxy = 0.0;
				double sumx = 0.0;
				double sumy = 0.0;
				double sumxsq = 0.0;
				Iterator<Long> workIter = workLog.iterator();
				Iterator<Long> timeIter = timeMsLog.iterator();
				while (workIter.hasNext()) {
					long work = workIter.next();
					long timeMs = timeIter.next();
					sumxy += work*timeMs;
					sumx += work;
					sumy += timeMs;
					sumxsq += work*work;
				}
				
				// solve for slope (a) and intercept (b)
				double a = (sumxy - sumx*sumy/workLog.size())/(sumxsq - sumx*sumx/workLog.size());
				double b = (sumy - a*sumx)/workLog.size();
				
				// extrapolate the finish time (y = ax+b), then compute the ETA
				double x = totalWork;
				return (long) (a*x + b) - timeMsLog.getLast();
			}
		},
		Quadratic {
			@Override
			public long getEtaMs(Deque<Long> workLog, Deque<Long> timeMsLog, long totalWork) {
				
				// NOTE: code shamelessly adapted from:
				// http://www.codeproject.com/KB/recipes/QuadraticRegression.aspx
				
				// not enough data?
				if (workLog.size() < 2) {
					return -1;
				}
				
				// compute our sums
				double s00 = workLog.size();
				double s10 = 0.0; // x
				double s20 = 0.0; // x^2
				double s30 = 0.0; // x^3
				double s40 = 0.0; // x^4
				double s01 = 0.0; // y
				double s11 = 0.0; // xy
				double s21 = 0.0; // x^2y
				Iterator<Long> workIter = workLog.iterator();
				Iterator<Long> timeIter = timeMsLog.iterator();
				while (workIter.hasNext()) {
					long work = workIter.next();
					long timeMs = timeIter.next();
					double x = work;
					double y = timeMs;
					s10 += x;
					s01 += y;
					s11 += x*y;
					x *= work;
					s20 += x;
					s21 += x*y;
					x *= work;
					s30 += x;
					x *= work;
					s40 += x;
				}
				
				// compute the quadratic model (y = ax^2 + bx + c)
				// UNDONE: if we really want to optimize this, we can pre-compute the multiplications
				double a = (
					s21 * ( s20 * s00 - s10 * s10 )
					- s11 * ( s30 * s00 - s10 * s20 )
					+ s01 * ( s30 * s10 - s20 * s20 )
				) / (
					s40 * ( s20 * s00 - s10 * s10 )
					- s30 * ( s30 * s00 - s10 * s20 )
					+ s20 * ( s30 * s10 - s20 * s20 )
				);
				
				double b = (
					s40 * ( s11 * s00 - s01 * s10 )
					- s30 * ( s21 * s00 - s01 * s20 )
					+ s20 * ( s21 * s10 - s11 * s20 )
				) / (
					s40 * ( s20 * s00 - s10 * s10 )
					- s30 * ( s30 * s00 - s10 * s20 )
					+ s20 * ( s30 * s10 - s20 * s20 )
				);
				
				double c = (
					s40 * ( s20 * s01 - s10 * s11 )
					- s30 * ( s30 * s01 - s10 * s21 )
					+ s20 * ( s30 * s11 - s20 * s21 )
				) / (
					s40 * ( s20 * s00 - s10 * s10 )
					- s30 * ( s30 * s00 - s10 * s20 )
					+ s20 * ( s30 * s10 - s20 * s20 )
				);
				
				// extrapolate the finish time, then compute the ETA
				double x = totalWork;
				return (long) (a*x*x + b*x + c) - timeMsLog.getLast();
			}
		};
		
		public abstract long getEtaMs(Deque<Long> workLog, Deque<Long> timeMsLog, long totalWork);
	}
	
	private long totalWork;
	private long reportIntervalMs;
	private Model model;
	private long currentWork;
	private Stopwatch stopwatch;
	private Deque<Long> workLog;
	private Deque<Long> timeMsLog;
	private long lastReportMs;
	private boolean reportMemory;
	
	public Progress(long totalWork) {
		this(totalWork, DefaultReportIntervalMs, DefaultModel);
	}
	
	public Progress(long totalWork, long reportIntervalMs) {
		this(totalWork, reportIntervalMs, DefaultModel);
	}
	
	public Progress(long totalWork, long reportIntervalMs, Model model) {
		
		// save params
		this.totalWork = totalWork;
		this.reportIntervalMs = reportIntervalMs;
		this.model = model;
		
		// init defaults
		currentWork = 0;
		stopwatch = new Stopwatch();
		workLog = new ArrayDeque<>();
		timeMsLog = new ArrayDeque<>();
		lastReportMs = 0;
		reportMemory = false;
		
		// add the 0,0 point to the work log
		stopwatch.start();
		addLog(0, 0);
	}
	
	public void setReportMemory(boolean val) {
		reportMemory = val;
	}
	
	public long getNumWorkDone() {
		return currentWork;
	}
	
	public long getTotalWork() {
		return totalWork;
	}
	public void setTotalWork(long val) {
		totalWork = val;
	}
	
	public boolean isFinished() {
		return currentWork == totalWork;
	}
	
	public boolean setProgress(long currentWork) {
		
		this.currentWork = currentWork;
		
		// should we update and/or report?
		long elapsedMs = stopwatch.getTimeNs()/NSpMS;
		boolean update = elapsedMs - timeMsLog.getLast() >= UpdateIntervalMs;
		boolean shouldReport = elapsedMs - lastReportMs >= reportIntervalMs;
		
		// if this is the last work done, force an update and a report
		if (isFinished()) {
			stopwatch.stop();
			update = true;
			shouldReport = true;
		}
		
		// update the work log if needed
		if (update) {
			addLog(currentWork, elapsedMs);
		}
		
		// report the progress if needed
		if (shouldReport) {
			
			// calculate the time remaining
			long etaMs = model.getEtaMs(workLog, timeMsLog, totalWork);
			
			// the extrapolation models aren't perfect, hide any errors
			if (etaMs < 0) {
				etaMs = 0;
			}
			
			System.out.print(String.format("Progress: %6.1f%%   ETA: %s",
				100.0*currentWork/totalWork,
				TimeFormatter.format(etaMs*NSpMS, 1)
			));
			if (reportMemory) {
				MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
				Runtime runtime = Runtime.getRuntime();
				System.out.print(String.format("   heapMem:%s, extMem:%s",
					JvmMem.getOldPool(),
					ExternalMemory.getUsageReport()
				));
			}
			System.out.println();
			lastReportMs = elapsedMs;
		}
		
		if (update) {
			pruneLog();
		}
		
		// add the finished message if needed
		if (currentWork == totalWork) {
			System.out.println("Finished in " + TimeFormatter.format(stopwatch.getTimeNs(), 1));
		}
		
		return shouldReport;
	}
	
	public boolean incrementProgress() {
		return incrementProgress(1);
	}
	
	public boolean incrementProgress(int numWorkDone) {
		return setProgress(currentWork + numWorkDone);
	}
	
	public long getTimeNs() {
		return stopwatch.getTimeNs();
	}
	
	private void addLog(long work, long timeMs) {
		workLog.add(work);
		timeMsLog.add(timeMs);
	}
	
	private void pruneLog() {
		
		// remove old entries from the log if needed
		while (workLog.size() > MaxLogSize) {
			workLog.removeFirst();
			timeMsLog.removeFirst();
		}
		assert (workLog.size() == timeMsLog.size());
	}
}
