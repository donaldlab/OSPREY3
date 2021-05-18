package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.AbstractQueuedSynchronizer;


/** a cyclic signal, ie it can be reused multiple times */
public class PulseSignal {

	private static class Sync extends AbstractQueuedSynchronizer {

		Sync() {
			setState(0);
		}

		@Override
		protected boolean tryAcquire(int acquires) {
			while (true) {
				int state = getState();
				if (state <= 0) {
					return false;
				}
				if (compareAndSetState(state, state - 1)) {
					return true;
				}
			}
		}

		@Override
		protected boolean tryRelease(int releases) {
			while (true) {
				int state = getState();
				if (compareAndSetState(state, state + 1)) {
					return true;
				}
			}
		}
	}

	private final Sync sync = new Sync();

	public void signal() {
		sync.release(1);
	}

	public enum Result {
		Signaled,
		TimedOut
	}

	public Result waitForSignal(long timeout, TimeUnit timeUnit) {
		try {
			if (sync.tryAcquireNanos(1, timeUnit.toNanos(timeout))) {
				return Result.Signaled;
			} else {
				return Result.TimedOut;
			}
		} catch (InterruptedException ex) {
			throw new Error(ex);
		}
	}
}
