package edu.duke.cs.osprey.coffee;

import edu.duke.cs.osprey.tools.TimeFormatter;

import java.time.Duration;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Consumer;


/**
 * Collects statistics on node processing.
 * Each instance of this class is designed to run on a separate thread.
 * All the threads communicated by eventually synchronizing to the Shared instance.
 */
public class NodeStats {

	private static final Object lock = new Object();

	public static class Values {

		public long rescored = 0;
		public long minimized = 0;
		public long expanded = 0;

		public void copyTo(Values other) {
			other.rescored = rescored;
			other.minimized = minimized;
			other.expanded = expanded;
		}

		public void addTo(Values other) {
			other.rescored += rescored;
			other.minimized += minimized;
			other.expanded += expanded;
		}

		public void clear() {
			rescored = 0;
			minimized = 0;
			expanded = 0;
		}
	}

	/**
	 * A single thread's view of the stats.
	 * Handles synchronization with other threads automatically.
	 */
	public class ForThread {

		private final Values localValues = new Values();

		private long lastSyncNs = 0;

		public void rescored() {
			localValues.rescored += 1;
			maybeSync();
		}

		public void minimized() {
			localValues.minimized += 1;
			maybeSync();
		}

		public void expanded() {
			localValues.expanded += 1;
			maybeSync();
		}

		private void maybeSync() {

			// if we've synced recently, don't bother
			long nowNs = System.nanoTime();
			long elapsedNs = nowNs - lastSyncNs;
			if (elapsedNs < syncIntervalNs) {
				return;
			}

			lastSyncNs = nowNs;
			sync();
		}

		public void sync() {
			synchronized (lock) {
				localValues.addTo(globalValues);
			}
			localValues.clear();
		}
	}

	public static class Report {

		public final long startNs;
		public final long stopNs;
		public final Values values;

		public Report(long startNs, long stopNs, Values values) {
			this.startNs = startNs;
			this.stopNs = stopNs;
			this.values = values;
		}

		@Override
		public String toString() {
			return String.format("NodeStats[res=%d, min=%d, exp=%d, in %s]",
				values.rescored,
				values.minimized,
				values.expanded,
				TimeFormatter.format(stopNs - startNs)
			);
		}
	}

	private long syncIntervalNs = Duration.ofSeconds(5).toNanos();

	private final Values globalValues = new Values();
	private long globalStartNs = 0;

	public void start() {
		globalStartNs = System.nanoTime();
	}

	public void start(Duration syncInterval) {
		start();
		syncIntervalNs = syncInterval.toNanos();
	}

	public Report get() {

		var values = new Values();
		long startNs;
		long stopNs;

		synchronized (lock) {
			startNs = globalStartNs;
			stopNs = System.nanoTime();
			globalStartNs = stopNs;
			globalValues.copyTo(values);
			globalValues.clear();
		}

		return new Report(startNs, stopNs, values);
	}

	/**
	 * Starts a thread to periodically report node stats
	 */
	public class Reporter implements Runnable, AutoCloseable {

		public final Duration interval;
		public final Consumer<Report> callback;

		private final AtomicBoolean isRunning = new AtomicBoolean(true);

		public Reporter(Duration interval, Consumer<Report> callback) {

			this.interval = interval;
			this.callback = callback;

			var thread = new Thread(this);
			thread.setName("NodeStatsReporter");
			thread.setDaemon(true);
			thread.start();
		}

		@Override
		public void run() {

			while (isRunning.get()) {

				// wait a bit
				try {
					// chill, linter... busy waiting is fine here
					//noinspection BusyWait
					Thread.sleep(interval.toMillis());
				} catch (InterruptedException ex) {
					break;
				}

				// do a report
				callback.accept(get());
			}
		}

		@Override
		public void close() {
			isRunning.set(false);
		}
	}
}
