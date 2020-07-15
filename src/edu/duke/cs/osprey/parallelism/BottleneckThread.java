package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Supplier;


/**
 * A worker thread that forces all work submitted from other threads to be serialized.
 *
 * Useful for when many threads need to access resources owned by only one thread.
 *
 * Provides shared access to thread-local stuff, thread-localotally prevents races, but not high-throughput.
 */
public class BottleneckThread implements AutoCloseable {

	public final String name;

	private final SynchronousQueue<Runnable> queue;
	private final Thread thread;

	private volatile boolean isRunning;

	public BottleneckThread(String name) {

		this.name = name;

		isRunning = true;

		queue = new SynchronousQueue<>();

		thread = new Thread(() -> threadLoop());
		thread.setName(name);
		thread.setDaemon(false);
		thread.start();
	}

	@Override
	public void close() {
		isRunning = false;
		try {
			thread.join();
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	private void threadLoop() {
		try {

			while (isRunning) {

				Runnable task = queue.poll(500, TimeUnit.MILLISECONDS);
				if (task != null) {
					task.run();
				}
			}

		} catch (InterruptedException ex) {
			ex.printStackTrace(System.err);
			// just exit the thread
		}
	}

	/**
	 * Run the task on the thread, and don't wait for it to finish.
	 */
	public void launch(Runnable task) {

		// don't allow tasks from the bottleneck thread, it'll cause a deadlock
		if (Thread.currentThread() == thread) {
			throw new IllegalArgumentException("can't launch new tasks from the bottleneck thread itself!");
		}

		try {
			while (true) {
				if (!isRunning) {
					throw new IllegalStateException("BottleneckThread is not running");
				}
				boolean wasAdded = queue.offer(task, 500, TimeUnit.MILLISECONDS);
				if (wasAdded) {
					break;
				}
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	/**
	 * Run the task on the thread, and wait for it to finish.
	 */
	public void exec(Runnable task) {

		var latch = new CountDownLatch(1);

		launch(() -> {
			try {
				task.run();
			} finally {
				latch.countDown();
			}
		});

		// wait for the task to finish
		try {
			while (true) {
				boolean finished = latch.await(500, TimeUnit.MILLISECONDS);
				if (finished) {
					break;
				}
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}
	}

	/**
	 * Run the task on the thread, wait for it to finish, and return the result.
	 */
	public <T> T get(Supplier<T> task) {

		var ref = new AtomicReference<T>(null);

		exec(() -> ref.set(task.get()));

		return ref.get();
	}
}
