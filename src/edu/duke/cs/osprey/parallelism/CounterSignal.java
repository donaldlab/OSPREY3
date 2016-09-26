package edu.duke.cs.osprey.parallelism;

public class CounterSignal {
	
	public static interface SignalCondition {
		boolean shouldSignal(int counter);
	}
	
	private int count;
	private SignalCondition condition;
	
	public CounterSignal(int initialCount, SignalCondition condition) {
		this.count = initialCount;
		this.condition = condition;
	}
	
	public synchronized void waitForSignal() {
		waitForSignal(0);
	}
	
	public synchronized void waitForSignal(long timeoutMs) {
		if (!condition.shouldSignal(count)) {
			try {
				wait(timeoutMs);
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
	
	public synchronized void offset(int delta) {
		count += delta;
		if (condition.shouldSignal(count)) {
			notifyAll();
		}
	}
}
