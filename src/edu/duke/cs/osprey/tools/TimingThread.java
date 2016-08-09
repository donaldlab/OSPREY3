package edu.duke.cs.osprey.tools;

import java.util.List;

public abstract class TimingThread extends Thread {
	
	private boolean isWarmedUp;
	private Object warmupSignal;
	private Object goSignal;
	
	public TimingThread(String name) {
		super(name);
		setDaemon(true);
		
		isWarmedUp = false;
		warmupSignal = new Object();
		goSignal = new Object();
	}
	
	protected abstract void warmup();
	protected abstract void time();
	
	@Override
	public void run() {
		
		warmup();
		isWarmedUp = true;
		
		synchronized (warmupSignal) {
			warmupSignal.notify();
		}
		
		synchronized (goSignal) {
			try {
				goSignal.wait();
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
		
		time();
	}
	
	public void waitForWarmup() {
		synchronized (warmupSignal) {
			if (!isWarmedUp) {
				try {
					warmupSignal.wait();
				} catch (InterruptedException ex) {
					throw new Error(ex);
				}
			}
		}
	}
	
	public void go() {
		synchronized (goSignal) {
			goSignal.notify();
		}
	}
	
	public static Stopwatch timeThreads(List<TimingThread> threads) {
	
		// start the threads and do the warmup
		for (TimingThread thread : threads) {
			thread.start();
		}
		
		// wait for the warmup to finish
		for (TimingThread thread : threads) {
			thread.waitForWarmup();
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// start the main event
		for (TimingThread thread : threads) {
			thread.go();
		}
		
		// wait for it to finish
		for (TimingThread thread : threads) {
			try {
				thread.join();
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
		
		stopwatch.stop();
		return stopwatch;
	}
}
