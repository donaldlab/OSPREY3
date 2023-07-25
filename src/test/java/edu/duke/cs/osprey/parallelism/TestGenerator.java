package edu.duke.cs.osprey.parallelism;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.concurrent.CountDownLatch;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class TestGenerator {

	@Test
	public void countTo3() {

		try (var gen = new Generator<Integer>(yielder -> {
			for (int i=0; i<3; i++) {
				yielder.yield(i);
			}
		})) {
			var iter = gen.iterator();
			assertThat(iter.hasNext(), is(true));
			assertThat(iter.next(), is(0));
			assertThat(iter.hasNext(), is(true));
			assertThat(iter.next(), is(1));
			assertThat(iter.hasNext(), is(true));
			assertThat(iter.next(), is(2));
			assertThat(iter.hasNext(), is(false));
		}
	}

	@Test
	public void iterateTooFar() {
		Assertions.assertThrows(NoSuchElementException.class, () -> {
			try (var gen = new Generator<Integer>(yielder -> {
				for (int i=0; i<3; i++) {
					yielder.yield(i);
				}
			})) {
				var iter = gen.iterator();
				iter.next();
				iter.next();
				iter.next();
				iter.next();
			}
		});
	}

	@Test
	public void threadFail() {

		Assertions.assertThrows(Generator.GeneratorFailedException.class, () -> {
			try (var gen = new Generator<Integer>(yielder -> {
				throw new Error("fail");
			})) {
				var iter = gen.iterator();
				iter.next();
			}
		});
	}

	@Test
	public void countTo100() {

		try (var gen = new Generator<Integer>(yielder -> {
			for (int i=0; i<100; i++) {
				yielder.yield(i);
			}
		})) {
			var iter = gen.iterator();
			for (int i=0; i<100; i++) {
				assertThat(iter.hasNext(), is(true));
				assertThat(iter.next(), is(i));
			}
			assertThat(iter.hasNext(), is(false));
		}
	}

	@Test
	public void countTo100_2threads() {

		final int numInts = 100;
		final int numThreads = 2;

		try (var gen = new Generator<Integer>(yielder -> {
			for (int i=0; i<numInts; i++) {
				yielder.yield(i);
			}
		})) {
			var iter = gen.iterator();

			var ints = new HashSet<Integer>();

			var latch = new CountDownLatch(1);

			// read from multiple threads
			var threads = IntStream.range(0, numThreads)
				.mapToObj(threadi -> new Thread(() -> {

					// wait for all the threads to be ready
					try {
						latch.await();
					} catch (InterruptedException ex) {
						throw new Error(ex);
					}

					// race on all the reads
					for (int i : gen) {
						ints.add(i);
					}
				}))
				.collect(Collectors.toList());

			// run the threads and wait for them to finish
			for (var thread : threads) {
				thread.start();
			}
			latch.countDown();
			try {
				for (var thread : threads) {
					thread.join();
				}
			} catch (InterruptedException ex) {
				throw new Error(ex);
			}

			// finally, we should have counted each int exactly once
			assertThat(ints.size(), is(numInts));
		}
	}
}
