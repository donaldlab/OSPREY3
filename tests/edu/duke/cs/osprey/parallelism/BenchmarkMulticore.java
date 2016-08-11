package edu.duke.cs.osprey.parallelism;

import edu.duke.cs.osprey.tools.Stopwatch;

public class BenchmarkMulticore {
	
	public static void main(String[] args)
	throws Exception {
		
		// time on the main thread
		System.out.println("Benchmarking main thread...");
		Stopwatch stopwatch = new Stopwatch().start();
		doWork();
		stopwatch.stop();
		System.out.println("Finished in " + stopwatch.getTime(2));
		long baseNs = stopwatch.getTimeNs();
		
		int[] numThreadsList = { 1, 2, 4, 8, 16, 32 };
		for (int numThreads : numThreadsList) {
			go(numThreads, baseNs);
		}
	}
	
	private static void doWork() {
		
		// just spin, it's the simplest possible workload...
		for (long i=0; i<1e9; i++);
	}
	
	private static void doWorkWithTiming(long baseNs) {
		Stopwatch stopwatch = new Stopwatch().start();
		
		doWork();
		
		stopwatch.stop();
		System.out.println(String.format("\t%s finished in %s, speed: %.1f%%",
			Thread.currentThread().getName(),
			stopwatch.getTime(2),
			100f*baseNs/stopwatch.getTimeNs()
		));
	}
	
	private static void go(int numThreads, long baseNs)
	throws Exception {
		
		System.out.println("Benchmarking " + numThreads + " threads...");
		
		// make the simplest possible multi-threaded work system
		
		Thread[] threads = new Thread[numThreads];
		for (int i=0; i<numThreads; i++) {
			threads[i] = new Thread("Spin-" + i) {
				@Override
				public void run() {
					doWorkWithTiming(baseNs);
				}
			};
		}
		
		Stopwatch stopwatch = new Stopwatch().start();
		
		for (Thread thread : threads) {
			thread.start();
		}
		for (Thread thread : threads) {
			thread.join();
		}
		
		/* try java 8's parallel streams
		List<Integer> numbers = new ArrayList<>();
		for (int i=0; i<numThreads; i++) {
			numbers.add(i);
		}
		numbers.parallelStream().forEach(i -> doWorkWithTiming(baseNs));
		*/
		
		stopwatch.stop();
		System.out.println(String.format("Finished in %s, speedup: %.2fx",
			stopwatch.getTime(2),
			(float)baseNs*numThreads/stopwatch.getTimeNs()
		));
	}
}
