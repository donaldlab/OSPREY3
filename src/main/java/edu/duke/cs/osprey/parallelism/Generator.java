package edu.duke.cs.osprey.parallelism;

import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.function.Consumer;


/**
 * Uses a thread to emulate a continuation-powered generator function.
 * Always use with try-with-resources to prevent leaking the associated thread!
 *
 * Only one iterator will ever be created for each generator and that iterator is thread-safe.
 */
public class Generator<T> implements Iterable<T>, AutoCloseable {

	public interface Yielder<T> {
		void yield(T value);
	}

	private final Thread genThread;
	private final Iterator<T> iter;

	private final PulseSignal iterSignal = new PulseSignal();
	private final SynchronousQueue<T> queue = new SynchronousQueue<>();
	private final AtomicBoolean isActive = new AtomicBoolean(true);

	/** a unique value to signal the end of the generator */
	@SuppressWarnings("unchecked")
	private final T end = (T)new Object();

	private static class AbortGeneratorException extends RuntimeException {}

	public static class GeneratorFailedException extends RuntimeException {
		public GeneratorFailedException(String msg) {
			super(msg);
		}
	}

	public Generator(Consumer<Yielder<T>> func) {

		// start the generator thread
		genThread = new Thread(() -> {

			// wait for the first iterator signal
			while (iterSignal.waitForSignal(1, TimeUnit.SECONDS) == PulseSignal.Result.TimedOut) {
				if (!isActive.get()) {
					return;
				}
			}

			var yielder = new Yielder<T>() {
				@Override
				public void yield(T value) {
					try {

						// send back the generated value
						while (!queue.offer(value, 1, TimeUnit.SECONDS)) {
							if (!isActive.get()) {
								throw new AbortGeneratorException();
							}
						}

						// wait for the next signal from the iterator
						while (iterSignal.waitForSignal(1, TimeUnit.SECONDS) == PulseSignal.Result.TimedOut) {
							if (!isActive.get()) {
								return;
							}
						}

					} catch (InterruptedException ex) {
						throw new Error(ex);
					}
				}
			};

			// call the generator function
			try {
				func.accept(yielder);

				// generator finished, send the end value
				while (!queue.offer(end, 1, TimeUnit.SECONDS)) {
					if (!isActive.get()) {
						throw new AbortGeneratorException();
					}
				}

			} catch (AbortGeneratorException ex) {
				// exit the thread
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}
		});
		genThread.setName("Generator Thread");
		genThread.setDaemon(true);
		genThread.start();

		// create the iterator, only one is allowed for this generator
		iter = new Iterator<T>() {

			private T nextValue = get();

			T get() {

				// tell the generator function to make the next value
				iterSignal.signal();

				// wait for it to show up
				try {
					while (true) {
						T value = queue.poll(1, TimeUnit.SECONDS);
						if (value != null) {

							// look for the special end value
							if (value == end) {
								return null;
							}

							// otherwise, return the generated value
							return value;
						}

						// timed out, make sure the generator thread is still active
						if (!genThread.isAlive()) {
							throw new GeneratorFailedException("Generator thread failed, no value can be iterated");
						}
					}
				} catch (InterruptedException ex) {
					throw new Error(ex);
				}
			}

			@Override
			public synchronized boolean hasNext() {
				return nextValue != null;
			}

			@Override
			public synchronized T next() {
				T value = nextValue;
				if (value == null) {
					throw new NoSuchElementException();
				}
				nextValue = get();
				return value;
			}
		};
	}

	@Override
	public void close() {
		isActive.set(false);
	}

	@Override
	public Iterator<T> iterator() {
		return iter;
	}
}
