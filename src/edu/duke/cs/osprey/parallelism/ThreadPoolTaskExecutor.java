package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import edu.duke.cs.osprey.tools.Cleaner;
import edu.duke.cs.osprey.tools.Cleaner.Cleanable;
import edu.duke.cs.osprey.tools.Cleaner.GarbageDetectable;

public class ThreadPoolTaskExecutor extends TaskExecutor implements GarbageDetectable {
	
	private static final ThreadFactory DaemonThreadFactory = new ThreadFactory() {
		// here there be daemons!
		
		private ThreadFactory threadFactory = Executors.defaultThreadFactory();

		@Override
		public Thread newThread(Runnable runnable) {
			Thread thread = threadFactory.newThread(runnable);
			thread.setDaemon(true);
			return thread;
		}
	};
	
	private static class Threads implements Cleanable {
		
		ThreadPoolExecutor pool;
		ThreadPoolExecutor listener;
		BlockingQueue<Runnable> queue;
		
		public Threads(int numThreads, int queueSize) {
			
			if (queueSize <= 0) {
				queue = new SynchronousQueue<>();
			} else {
				queue = new ArrayBlockingQueue<>(queueSize);
			}
			pool = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.DAYS, queue, DaemonThreadFactory);
			pool.prestartAllCoreThreads();
			
			// use an unbounded queue for the listener thread
			// let task results pile up until the listener thread can process them
			listener = new ThreadPoolExecutor(1, 1, 0, TimeUnit.DAYS, new LinkedBlockingQueue<>(), DaemonThreadFactory);
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
	}
	
	/**
	 * Controls task queue size.
	 * Set this to 0 to cause main thread to block until task thread is ready.
	 * Set this to >0 to "buffer" tasks so task threads don't have to wait on the main thread to start a task.
	 * >0 can sometimes be faster than 0, but only works if you know how many tasks you have in advance.
	 * otherwise, you can end up execting more tasks than you need.
	 * The best queue size to use is determined by the amount of work it takes to create a task vs execute it.
	 * Experiment to find the best values for your problem.
	 */
	public int queueSize = 0;
	
	private Threads threads;
	private AtomicLong numTasksStarted;
	private AtomicLong numTasksFinished;
	private Signal taskSignal;
	
	public ThreadPoolTaskExecutor() {
		threads = null;
		numTasksStarted = new AtomicLong(0);
		numTasksFinished = new AtomicLong(0);
		taskSignal = new Signal();
	}
	
	public void start(int numThreads) {
		threads = new Threads(numThreads, queueSize);
		Cleaner.addCleaner(this, threads);
	}
	
	public void stop() {
		if (threads != null) {
			threads.clean();
			threads = null;
		}
	}
	
	public void stopAndWait(int timeoutMs) {
		if (threads != null) {
			threads.cleanAndWait(timeoutMs);
			threads = null;
		}
	}
	
	@Override
	public int getParallelism() {
		return threads.pool.getCorePoolSize();
	}
	
	@Override
	public void submit(Runnable task, TaskListener listener) {
		try {
			
			boolean wasAdded = false;
			while (!wasAdded) {
				
				// NOTE: don't use ThreadPoolExecutor.submit() to send tasks, because it won't let us block.
				// access the work queue directly instead, so we can block if the thread pool isn't ready yet.
				wasAdded = threads.queue.offer(() -> {
					
					// run the task
					task.run();
					
					// send the result to the listener thread
					threads.listener.submit(() -> {
						
						// run the listener
						listener.onFinished(task);
						
						// tell anyone waiting that we finished a task
						numTasksFinished.incrementAndGet();
						taskSignal.sendSignal();
					});
					
				}, 1, TimeUnit.SECONDS);
			}
			
			numTasksStarted.incrementAndGet();
			
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}
	
	@Override
	public void waitForFinish() {
		
		long numTasks = numTasksStarted.get();
		
		while (numTasksFinished.get() < numTasks) {
			
			// wait a bit before checking again, unless a task finishes
			taskSignal.waitForSignal(100);
		}
	}
}
