package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.TimeUnit;


/**
 * Runs a task no more than once every interval.
 */
public class RateLimitedThread implements AutoCloseable {

	public final String name;
	public final long interval;
	public final TimeUnit unit;
	public final Runnable task;

	private final Thread thread;

	private volatile boolean isRunning;
	private volatile boolean isRequested;

	public RateLimitedThread(String name, long interval, TimeUnit unit, Runnable task) {

		this.name = name;
		this.interval = interval;
		this.unit = unit;
		this.task = task;

		isRunning = true;
		isRequested = false;

		thread = new Thread(() -> threadLoop());
		thread.setName(name);
		thread.setDaemon(false);
		thread.start();
	}

	@Override
	public void close() {
		isRunning = false;
		try {
			thread.join();
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	public void request() {
		isRequested = true;
	}

	@SuppressWarnings("BusyWait")
	/*
		IntelliJ's linter doesn't like the sleep() call in the loop.
		Usually, it would be right. Polling is inefficient if you're waiting for an even to happen.
		But in this case, the whole point is to wait a specific amount of time, so this usage is fine.
	*/
	private void threadLoop() {

		long millis = unit.toMillis(interval);

		try {

			while (isRunning) {
				if (isRequested) {
					task.run();
					isRequested = false;
				}
				Thread.sleep(millis);
			}

		} catch (InterruptedException ex) {
			ex.printStackTrace(System.err);
			// just exit the thread
		}
	}
}
