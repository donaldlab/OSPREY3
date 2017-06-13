package edu.duke.cs.osprey.parallelism;

import edu.duke.cs.tpie.Cleaner.Cleanable;

public abstract class WorkThread extends Thread implements Cleanable {
	
	// this flag is hit from multiple threads concurrently, so make it volatile
	private volatile boolean isRunning;
	
	protected WorkThread(String name) {
		super(name);
		setDaemon(true);
		isRunning = false;
	}
	
	public boolean isRunning() {
		return isRunning();
	}
	
	public void askToStop() {
		isRunning = false;
	}
	
	public void askToStopAndWait()
	throws InterruptedException {
		askToStop();
		join();
	}
	
	@Override
	public void run() {
		
		isRunning = true;
		
		while (isRunning) {
			try {
				doWork();
			} catch (InterruptedException ex) {
				break;
			}
		}
		
		isRunning = false;
	}
	
	protected abstract void doWork() throws InterruptedException;
	
	@Override
	public void clean() {
		askToStop();
	}
}
