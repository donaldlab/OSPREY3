package edu.duke.cs.osprey.parallelism;

public class Signal {
	
	private boolean isSignaled = false;
	
	public synchronized void waitForSignal() {
		waitForSignal(0);
	}
	
	public synchronized void waitForSignal(long timeoutMs) {
		if (!isSignaled) {
			try {
				wait(timeoutMs);
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
	
	public synchronized void sendSignal() {
		isSignaled = true;
		notifyAll();
	}
	
	public synchronized void reset() {
		isSignaled = false;
	}
}
