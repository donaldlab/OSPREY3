package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicLong;


/**
 * Allows threads to wait for a large amount of work to be finished.
 */
public class WorkLatch {

	private final AtomicLong incoming;
	private final CountDownLatch latch;

	public WorkLatch(long size) {
		incoming = new AtomicLong(size);
		latch = new CountDownLatch(1);
	}

	public void await(long timeout, TimeUnit unit) {

		// anything to wait for?
		if (incoming.get() > 0) {

			// yup
			try {
				boolean finished = latch.await(timeout, unit);
				if (!finished) {
					// technically, the last work could come in *right after* we timed out,
					// but before this message gets constructed.
					// let's hope that doesn't happen much,
					// because it would be hard to fix with synchronization
					throw new RuntimeException("Timeed out waiting work to finish, remaining: " + incoming.get());
				}
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}
		}
	}

	public void finished(long size) {
		long remaining = incoming.updateAndGet(i -> i -= size);
		if (remaining <= 0) {
			latch.countDown();
			if (remaining < 0) {
				System.err.println("WARNING: too much work done, overwork: " + -remaining);
				Thread.dumpStack();
			}
		}
	}
}
