package edu.duke.cs.osprey.parallelism;

import java.util.List;

import edu.duke.cs.osprey.tools.Stopwatch;

public abstract class TimingThread extends Thread {
	
	public Signal warmupSignal;
	public Signal goSignal;
	
	public TimingThread(String name) {
		super(name);
		setDaemon(true);
		
		warmupSignal = new Signal();
		goSignal = new Signal();
	}
	
	protected abstract void warmup();
	protected abstract void time();
	
	@Override
	public void run() {
		warmup();
		warmupSignal.sendSignal();
		goSignal.waitForSignal();
		time();
	}
	
	public void waitForFinish() {
		try {
			join();
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}
	
	public static Stopwatch timeThreads(List<TimingThread> threads) {
	
		// start the threads and do the warmup
		for (TimingThread thread : threads) {
			thread.start();
		}
		
		// wait for the warmup to finish
		for (TimingThread thread : threads) {
			thread.warmupSignal.waitForSignal();
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		// start the main event
		for (TimingThread thread : threads) {
			thread.goSignal.sendSignal();
		}
		
		// wait for it to finish
		for (TimingThread thread : threads) {
			thread.waitForFinish();
		}
		
		stopwatch.stop();
		return stopwatch;
	}
}
