package edu.duke.cs.osprey.parallelism;

public class Signal {
	
	public boolean isSignaled = false;
	
	public synchronized void waitForSignal() {
		if (!isSignaled) {
			try {
				wait();
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
	
	public synchronized void sendSignal() {
		isSignaled = true;
		notifyAll();
	}
}
