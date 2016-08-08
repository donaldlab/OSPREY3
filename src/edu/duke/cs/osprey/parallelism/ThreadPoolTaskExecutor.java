package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

public class ThreadPoolTaskExecutor extends TaskExecutor {
	
	private static class ListeningTask implements Runnable {
		
		public Runnable task;
		public TaskExecutor.TaskListener listener;
		
		public ListeningTask(Runnable task, TaskExecutor.TaskListener listener) {
			this.task = task;
			this.listener = listener;
		}
		
		@Override
		public void run() {
			task.run();
		}
	}
	
	private BlockingQueue<Runnable> incomingQueue;
	private ThreadPoolExecutor pool;
	private Thread listenerThread;
	private BlockingQueue<Runnable> outgoingQueue;
	private int numTasksToFinish;
	private int numTasksFinished;
	private boolean isFinishing; // NOTE: we could make this volatile or not, it doesn't matter much
	private Object finishSignal;
	
	public ThreadPoolTaskExecutor() {
		incomingQueue = null;
		pool = null;
		listenerThread = null;
		outgoingQueue = null;
		numTasksToFinish = 0;
		numTasksFinished = 0;
		isFinishing = false;
		finishSignal = new Object();
	}
	
	public void start(int numThreads) {
		
		// make the incoming queue about twice the number of threads,
		// so the main thread can keep it as full as possible
		incomingQueue = new ArrayBlockingQueue<>(numThreads*2);
		
		pool = new ThreadPoolExecutor(numThreads, numThreads, 0, TimeUnit.SECONDS, incomingQueue) {
			
			@Override
			public void afterExecute(Runnable task, Throwable error) {
					
				// was there an error?
				if (error != null) {

					// update the finish detection
					synchronized (this) {
						numTasksFinished++;
					}
					
					// NOTE: the error already gets written to the console by the thread pool
					
				} else {
					
					try {
						
						// pass off to listener thread, wait if needed
						boolean wasAdded = false;
						while (!wasAdded) {
							wasAdded = outgoingQueue.offer(task, 1, TimeUnit.SECONDS);
						}
						
					} catch (InterruptedException ex) {
						throw new RuntimeException(ex);
					}
				}
			}
		};
		pool.setThreadFactory(new ThreadFactory() {
			@Override
			public Thread newThread(Runnable task) {
				
				// set the daemon status on the thread pool threads
				// so the thread pool shuts down automatically if all the main threads exit
				Thread thread = Executors.defaultThreadFactory().newThread(task);
				thread.setDaemon(true);
				return thread;
			}
		});
		pool.prestartAllCoreThreads();
		
		// init finish detection state
		numTasksToFinish = 0;
		numTasksFinished = 0;
		isFinishing = false;
		
		// start the listener thread if needed
		outgoingQueue = new ArrayBlockingQueue<>(numThreads*2);
		listenerThread = new Thread("TaskExecutor-listener") {
			
			// this is like a constructor =)
			{
				setDaemon(true);
			}
			
			@Override
			public void run() {
				try {
					while (!pool.isTerminated()) {
					
						// wait for the next task to finish, or timeout
						Runnable task = outgoingQueue.poll(1, TimeUnit.SECONDS);
						if (task != null) {
							
							// call the listener if needed
							if (task instanceof ListeningTask) {
								ListeningTask listeningTask = (ListeningTask)task;
								listeningTask.listener.onFinished(listeningTask.task);
							}
							
							// is this the last task?
							boolean isLastTask;
							synchronized (this) {
								numTasksFinished++;
								isLastTask = isFinishing && numTasksFinished == numTasksToFinish;
							}
							
							if (isLastTask) {
								
								// we're done!
								synchronized (this) {
									isFinishing = false;
								}
								
								// tell anyone who cares that we're done processing tasks
								synchronized (finishSignal) {
									finishSignal.notifyAll();
								}
							}
						}
					}
						
				} catch (InterruptedException ex) {
					// something told us to stop, so exit this thread
				}
			}
		};
		listenerThread.start();
	}
	
	@Override
	public void submit(Runnable task) {
		
		if (incomingQueue == null) {
			throw new IllegalStateException("thread pool is not started, can't submit tasks");
		}
		
		// update finish detection
		synchronized (this) {
			numTasksToFinish++;
		}
		
		// send the task, and wait for space in the queue if needed
		try {
			boolean wasAdded = false;
			while (!wasAdded) {
				wasAdded = incomingQueue.offer(task, 1, TimeUnit.SECONDS);
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	@Override
	public void submit(Runnable task, TaskListener listener) {
		submit(new ListeningTask(task, listener));
	}
	
	@Override
	public void waitForFinish() {
		try {
			isFinishing = true;
			synchronized (finishSignal) {
				finishSignal.wait();
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}
	
	public boolean stopAndWait(int timeoutMs)
	throws InterruptedException {
		pool.shutdown();
		if (!pool.awaitTermination(timeoutMs, TimeUnit.MILLISECONDS)) {
			return false;
		}
		listenerThread.join(timeoutMs);
		return outgoingQueue.isEmpty();
	}
}
