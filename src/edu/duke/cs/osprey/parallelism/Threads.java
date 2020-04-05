package edu.duke.cs.osprey.parallelism;

import edu.duke.cs.osprey.tools.AutoCleanable;

import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;


public class Threads implements AutoCleanable {

	private static int nextId = 0;

	private final int poolId = nextId++;
	private final ThreadPoolExecutor pool;
	private final ThreadPoolExecutor listener;

	public final BlockingQueue<Runnable> queue;

	public Threads(int numThreads, int queueSize) {

		if (queueSize <= 0) {
			queue = new SynchronousQueue<>();
		} else {
			queue = new ArrayBlockingQueue<>(queueSize);
		}
		AtomicInteger threadId = new AtomicInteger(0);
		pool = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.DAYS, queue, (runnable) -> {
			Thread thread = Executors.defaultThreadFactory().newThread(runnable);
			thread.setDaemon(true);
			thread.setName(String.format("pool-%d-%d", poolId, threadId.getAndIncrement()));
			return thread;
		});
		pool.prestartAllCoreThreads();

		// use an unbounded queue for the listener thread
		// let task results pile up until the listener thread can process them
		listener = new ThreadPoolExecutor(1, 1, 0, TimeUnit.DAYS, new LinkedBlockingQueue<>(), (runnable) -> {
			Thread thread = Executors.defaultThreadFactory().newThread(runnable);
			thread.setDaemon(true);
			thread.setName(String.format("pool-%d-listener", poolId));
			return thread;
		});
		listener.prestartAllCoreThreads();
	}

	@Override
	public void clean() {
		pool.shutdown();
		listener.shutdown();
	}

	public void cleanAndWait(int timeoutMs) {
		clean();
		try {
			pool.awaitTermination(timeoutMs, TimeUnit.MILLISECONDS);
			listener.awaitTermination(timeoutMs, TimeUnit.MILLISECONDS);
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}

	public int size() {
		return pool.getCorePoolSize();
	}

	public boolean submit(long timeout, TimeUnit timeUnit, Runnable runnable) {
		// NOTE: don't use ThreadPoolExecutor.submit() to send tasks, because it won't let us block.
		// access the work queue directly instead, so we can block if the thread pool isn't ready yet.
		try {
			return queue.offer(runnable, timeout, timeUnit);
		} catch (InterruptedException ex) {
			throw new RuntimeException("can't submit runnable to threads", ex);
		}
	}

	public void submitLoop(long interval, TimeUnit timeUnit, Runnable runnable) {
		boolean wasAdded = false;
		while (!wasAdded) {
			wasAdded = submit(interval, timeUnit, runnable);
		}
	}

	public void submitToListener(Runnable runnable) {
		listener.submit(runnable);
	}
}
