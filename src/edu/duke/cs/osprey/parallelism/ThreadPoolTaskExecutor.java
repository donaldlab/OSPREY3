package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.RejectedExecutionException;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;

import edu.duke.cs.osprey.tools.Cleaner;
import edu.duke.cs.osprey.tools.Cleaner.GarbageDetectable;

public class ThreadPoolTaskExecutor extends TaskExecutor implements GarbageDetectable {
	
	// TODO: use cleaner to call shutdown()
	// the JRE won't do it automatically when the threads are still running
	
	private ThreadPoolExecutor pool;
	private ThreadPoolExecutor listenerThread;
	
	private AtomicLong numTasksStarted;
	private AtomicLong numTasksFinished;
	
	public ThreadPoolTaskExecutor() {
		pool = null;
	}
	
	public void start(int numThreads) {
		
		// don't use queues for the main thread pool
		// handoff a task when a thread is ready
		pool = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.DAYS, new SynchronousQueue<>());
		pool.prestartAllCoreThreads();
		
		// use an unbounded queue for the listener thread
		// let task results pile up until the listener thread can process them
		listenerThread = new ThreadPoolExecutor(1, 1, 0, TimeUnit.DAYS, new LinkedBlockingQueue<>());
		listenerThread.prestartAllCoreThreads();
	
		Cleaner.addCleaner(this, () -> {
			if (pool != null) {
				pool.shutdown();
			}
			if (listenerThread != null) {
				listenerThread.shutdown();
			}
		});
		
		numTasksStarted = new AtomicLong(0);
		numTasksFinished = new AtomicLong(0);
	}
	
	public void stop() {
		pool.shutdown();
	}
	
	public void stopAndWait(int timeoutMs) {
		stop();
		try {
			pool.awaitTermination(timeoutMs, TimeUnit.MILLISECONDS);
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}
	
	@Override
	public int getParallelism() {
		return pool.getCorePoolSize();
	}
	
	@Override
	public void submit(Runnable task, TaskListener listener) {
		
		while (true) {
			try {
				
				pool.submit(() -> {
					
					// run the task
					task.run();
					
					// send the result to the listener thread
					listenerThread.submit(() -> {
						
						// run the listener
						listener.onFinished(task);
						
						// tell anyone waiting that we finished a task
						numTasksFinished.incrementAndGet();
						synchronized (numTasksFinished) {
							numTasksFinished.notifyAll();
						}
						
					});
				});
				
				numTasksStarted.incrementAndGet();
				
				// task submitted, we're done here
				break;
				
			} catch (RejectedExecutionException ex) {
				// keep trying to submit
			}
			
			// wait a bit before trying again
			try {
				Thread.sleep(20);
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
	
	@Override
	public void waitForFinish() {
		
		long numTasks = numTasksStarted.get();
		
		while (numTasksFinished.get() < numTasks) {
			
			// wait for a task to finish before checking again
			try {
				synchronized (numTasksFinished) {
					numTasksFinished.wait(20);
				}
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		}
	}
}
