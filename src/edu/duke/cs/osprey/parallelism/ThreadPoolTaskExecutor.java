package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;

import edu.duke.cs.tpie.Cleaner;
import edu.duke.cs.tpie.Cleaner.Cleanable;
import edu.duke.cs.tpie.Cleaner.GarbageDetectable;


public class ThreadPoolTaskExecutor extends TaskExecutor implements GarbageDetectable {
	
	private static class Threads implements Cleanable {
		
		private static int nextId = 0;
		
		final int poolId = nextId++;
		final ThreadPoolExecutor pool;
		final ThreadPoolExecutor listener;
		final BlockingQueue<Runnable> queue;
		
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
	}
	
	/**
	 * Controls task queue size.
	 * Set this to 0 to cause main thread to block until task thread is ready.
	 * Set this to >0 to "buffer" tasks so task threads don't have to wait on the main thread to start a task.
	 * >0 can sometimes be faster than 0, but only works if you know how many tasks you have in advance.
	 * otherwise, you can end up executing more tasks than you need.
	 * The best queue size to use is determined by the amount of work it takes to create a task vs execute it.
	 * Experiment to find the best values for your problem.
	 */
	public int queueSize = 0;
	
	private Threads threads;
	private AtomicLong numTasksStarted;
	private AtomicLong numTasksFinished;
	private AtomicReference<TaskException> exception;
	private Signal taskSignal;
	
	public ThreadPoolTaskExecutor() {
		threads = null;
		numTasksStarted = new AtomicLong(0);
		numTasksFinished = new AtomicLong(0);
		exception = new AtomicReference<>(null);
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
	public void clean() {
		stop();
	}
	
	@Override
	public int getParallelism() {
		return threads.pool.getCorePoolSize();
	}

	@Override
	public boolean isBusy() {
		return getNumRunningTasks() >= getParallelism();
	}
	
	@Override
	public boolean isWorking() {
		return getNumRunningTasks() > 0;
	}
	
	@Override
	public <T> void submit(Task<T> task, TaskListener<T> listener) {
		try {
			
			boolean wasAdded = false;
			while (!wasAdded) {
				
				// check for exceptions
				// NOTE: waitForFinish will throw the exception
				if (exception.get() != null) {
					waitForFinish();
				}
				
				// NOTE: don't use ThreadPoolExecutor.submit() to send tasks, because it won't let us block.
				// access the work queue directly instead, so we can block if the thread pool isn't ready yet.
				wasAdded = threads.queue.offer(() -> {
					
					try {
					
						// run the task
						T result = task.run();
						
						// send the result to the listener thread
						threads.listener.submit(() -> {
							
							try {
								
								// run the listener
								listener.onFinished(result);
								
							} catch (Throwable t) {
								recordException(task, listener, t);
							}
							
							// tell anyone waiting that we finished a task
							finishedTask();
						});
						
					} catch (Throwable t) {
						recordException(task, listener, t);
	
						// the task failed, but still report finish
						finishedTask();
					}
					
				}, 400, TimeUnit.MILLISECONDS);
			}
			
			// the task was started successfully, hooray!
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
		
		// check for exceptions
		TaskException t = exception.get();
		if (t != null) {
			throw t;
		}
	}
	
	public long getNumRunningTasks() {
		return numTasksStarted.get() - numTasksFinished.get();
	}
	
	private void recordException(Task<?> task, TaskListener<?> listener, Throwable t) {
		
		// record the exception, but don't overwrite any existing exceptions
		// TODO: keep a list of all exceptions?
		exception.compareAndSet(null, new TaskException(task, listener, t));
	}
	
	private void finishedTask() {
		numTasksFinished.incrementAndGet();
		taskSignal.sendSignal();
	}
}
