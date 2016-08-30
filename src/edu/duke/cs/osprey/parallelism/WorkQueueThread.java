package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;

public abstract class WorkQueueThread<T> extends WorkThread {
	
	public BlockingQueue<T> queue;
	
	public WorkQueueThread(String name, BlockingQueue<T> queue) {
		super(name);
		this.queue = queue;
	}
	
	@Override
	public void doWork()
	throws InterruptedException {
		
		// get the next piece of work in the queue
		T work = queue.poll(200, TimeUnit.MILLISECONDS);
		if (work != null) {
			doWork(work);
		}
	}
	
	protected abstract void doWork(T work) throws InterruptedException;
}
