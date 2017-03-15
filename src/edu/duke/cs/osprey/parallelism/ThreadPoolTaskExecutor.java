package edu.duke.cs.osprey.parallelism;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

import edu.duke.cs.osprey.tools.Cleaner;
import edu.duke.cs.osprey.tools.Cleaner.GarbageDetectable;

public class ThreadPoolTaskExecutor extends TaskExecutor implements TaskExecutor.NeedsCleanup, GarbageDetectable {
	
	private static class FailedTask implements Runnable {
		
		private Throwable error;
		
		public FailedTask(Throwable error) {
			this.error = error;
		}
		
		@Override
		public void run() {
			// don't need to do anything
		}
	}
	
	private static class FinishState extends Signal {
		
		private boolean submissionsOpen;
		private int numSubmitted;
		private int numWin;
		private int numFail;
		private Throwable error;
		
		public FinishState() {
			reset();
		}
		
		@Override
		public synchronized void reset() {
			super.reset();
			submissionsOpen = true;
			numSubmitted = 0;
			numWin = 0;
			numFail = 0;
			error = null;
		}
		
		public synchronized void recordSubmission() {
			if (!submissionsOpen) {
				reset();
			}
			numSubmitted++;
		}
		
		@Override
		public synchronized void waitForSignal() {
			submissionsOpen = false;
			
			// is there anything to wait for?
			if (numWin + numFail < numSubmitted) {
				super.waitForSignal();
			}
		}
	
		public synchronized void recordResult(Throwable error) {
			
			// record the result
			if (error != null) {
				numFail++;
				
				// forward the error to the main thread
				this.error = error;
				
			} else {
				numWin++;
			}
			
			// if this the last task, send the signal
			boolean isLastTask = numWin + numFail == numSubmitted;
			if (isLastTask) {
				sendSignal();
			}
		}
		
		public synchronized boolean hasError() {
			return error != null;
		}
		
		public synchronized void throwError() {
			if (error != null) {
				throw new RuntimeException("error while processing task", error);
			}
		}
	}
	
	private static class TaskThread extends WorkQueueThread<Runnable> {
		
		private int queueFactor;
		private CounterSignal idleThreadsSignal;
		private BlockingQueue<Runnable> outgoingQueue;
		
		private TaskThread(int index, BlockingQueue<Runnable> queue, BlockingQueue<Runnable> outgoingQueue, int queueFactor, CounterSignal idleThreadsSignal) {
			super("TaskThread-" + index, queue);
			this.queueFactor = queueFactor;
			this.idleThreadsSignal = idleThreadsSignal;
			this.outgoingQueue = outgoingQueue;
		}
		
		@Override
		public void doWork(Runnable task)
		throws InterruptedException {
			
			// update signals for waitForSpace()
			if (queueFactor == 0) {
				idleThreadsSignal.offset(-1);
			}
			
			// run the task
			try {
				task.run();
			} catch (Throwable t) {
				task = new FailedTask(t);
			}
			
			// signal task finished
			boolean wasAdded = false;
			while (!wasAdded) {
				wasAdded = outgoingQueue.offer(task, 1, TimeUnit.SECONDS);
			}
			
			// update signals for waitForSpace()
			if (queueFactor == 0) {
				idleThreadsSignal.offset(1);
			}
		}
	}
	
	private static class ListenerThread extends WorkQueueThread<Runnable> {
		
		private FinishState finishState;
		
		public ListenerThread(BlockingQueue<Runnable> queue, FinishState finishState) {
			super("TaskThread-listener", queue);
			this.finishState = finishState;
		}
	
		@Override
		public void doWork(Runnable task)
		throws InterruptedException {
			
			// was there a failure on a previous task?
			if (finishState.hasError()) {
				return;
			}
			
			// was there a failure on this task?
			Throwable error = null;
			if (task instanceof FailedTask) {
				error = ((FailedTask)task).error;
			}
			
			// should we call a listener?
			if (task instanceof ListeningTask) {
				try {
					
					ListeningTask listeningTask = (ListeningTask)task;
					listeningTask.listener.onFinished(listeningTask.task);
					
				} catch (Throwable t) {
					
					// the listener failed, report the error
					error = t;
				}
			}
			
			// record the result
			finishState.recordResult(error);
		}
	}
	
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
	
	private int queueFactor;
	private BlockingQueue<Runnable> incomingQueue;
	private BlockingQueue<Runnable> outgoingQueue;
	private List<TaskThread> threads;
	private ListenerThread listenerThread;
	private FinishState finishState;
	private CounterSignal idleThreadsSignal;
	
	public ThreadPoolTaskExecutor() {
		incomingQueue = null;
		outgoingQueue = null;
		threads = null;
		listenerThread = null;
		finishState = new FinishState();
		idleThreadsSignal = null;
	}
	
	public void start(int numThreads) {
		
		// by default, make the incoming queue equal to the number of threads,
		// so the worker threads probably never have to block waiting for the next task
		// the main thread will probably fill up the incoming queue with tasks while the task threads are running
		start(numThreads, 1);
		
		// NOTE: this works well for finite workloads (ie for loops)
		// for infinite workloads (ie in while loops), you'd probably want a queue factor of 0
		// to prevent waiting for extra tasks after the loop exits
	}
	
	public void start(int numThreads, int queueFactor) {
		
		// allocate queues
		this.queueFactor = queueFactor;
		incomingQueue = new ArrayBlockingQueue<>(Math.max(1, numThreads*queueFactor));
		outgoingQueue = new ArrayBlockingQueue<>(numThreads);
		
		// init state
		finishState.reset();
		
		// start all the threads
		threads = new ArrayList<>(numThreads);
		for (int i=0; i<numThreads; i++) {
			threads.add(Cleaner.addCleaner(this, new TaskThread(i, incomingQueue, outgoingQueue, queueFactor, idleThreadsSignal)));
		}
		for (TaskThread thread : threads) {
			thread.start();
		}
		
		// start the listener thread
		listenerThread = Cleaner.addCleaner(this, new ListenerThread(outgoingQueue, finishState));
		listenerThread.start();
		
		// init waitForSpace() signals
		idleThreadsSignal = new CounterSignal(numThreads, (int numIdleThreads) -> numIdleThreads > 0);
	}
	
	@Override
	public int getParallelism() {
		return threads.size();
	}
	
	@Override
	public void submit(Runnable task) {
		
		if (incomingQueue == null) {
			throw new IllegalStateException("thread pool is not started, can't submit tasks");
		}
		
		// send the task, and wait for space in the queue if needed
		try {
			finishState.recordSubmission(); // sync
			boolean wasAdded = false;
			while (!wasAdded) {
				wasAdded = incomingQueue.offer(task, 1, TimeUnit.SECONDS); // sync
				
				finishState.throwError();
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
	public void waitForSpace() {
		if (queueFactor == 0) {
			
			// wait for a task thread to become idle
			idleThreadsSignal.waitForSignal();
			
			finishState.throwError();
			
		} else {
			throw new Error("if you need to wait for space, you should set the queue factor to 0 in init()");
		}
	}
	
	@Override
	public void waitForFinish() {
		finishState.waitForSignal();
		finishState.throwError();
	}
	
	@Override
	public void cleanup() {
		stop();
	}
	
	public void stop() {
		askThreadsToStop();
		clearThreads();
	}
	
	public void stopAndWait(int timeoutMs)
	throws InterruptedException {
		askThreadsToStop();
		for (TaskThread thread : threads) {
			thread.join(timeoutMs);
		}
		listenerThread.join(timeoutMs);
		clearThreads();
	}
	
	private void askThreadsToStop() {
		for (TaskThread thread : threads) {
			thread.askToStop();
		}
		listenerThread.askToStop();
	}
	
	private void clearThreads() {
		threads.clear();
		listenerThread = null;
	}
}
