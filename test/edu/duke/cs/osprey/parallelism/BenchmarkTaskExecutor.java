package edu.duke.cs.osprey.parallelism;

import edu.duke.cs.osprey.tools.Stopwatch;

@SuppressWarnings("unused")
public class BenchmarkTaskExecutor {
	
	private static abstract class TaskFactory {
		
		public int numRuns;
		
		protected TaskFactory(int numRuns) {
			this.numRuns = numRuns;
		}
		
		public abstract Runnable makeTask();
	}
	
	public static void main(String[] args) {
	
		benchmarkBigCpu();
		//benchmarkManySmallCpu();
	}
	
	private static void benchmarkBigCpu() {
		
		class Task implements Runnable {
			@Override
			public void run() {
				
				// just spin
				for (int i=0; i<1e9; i++);
				
			}
		}
		
		benchmark(new TaskFactory(1) {
			@Override
			public Runnable makeTask() {
				return new Task();
			}
		});
	}
	
	private static void benchmarkManySmallCpu() {
		
		class Task implements Runnable {
			@Override
			public void run() {
				
				// just spin
				for (int i=0; i<1e6; i++);
				
			}
		}
		
		benchmark(new TaskFactory(2000) {
			@Override
			public Runnable makeTask() {
				return new Task();
			}
		});
	}
	
	private static void benchmark(TaskFactory factory) {
		
		int[] numThreadsList = { 1, 2, 4, 8, 16, 32 };
		ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
		
		long baseNs = Long.MAX_VALUE;
		
		for (int numThreads : numThreadsList) {
			for (int i=0; i<3; i++) {
		
				System.out.print(String.format("Benchmarking %2d threads...  ", numThreads));
				tasks.start(numThreads);
				
				// TEMP
				System.out.println();
				
				// do the warmup
				for (int j=0; j<numThreads; j++) {
					tasks.submit(factory.makeTask());
				}
				tasks.waitForFinish();
				
				// time the tasks
				Stopwatch stopwatch = new Stopwatch().start();
				for (int j=0; j<numThreads; j++) {
					tasks.submit(factory.makeTask());
				}
				tasks.waitForFinish();
				stopwatch.stop();
				
				// update timing info
				if (numThreads == 1) {
					baseNs = Math.min(baseNs, stopwatch.getTimeNs());
					System.out.println("Finished in " + stopwatch.getTime(2));
				} else {
					System.out.println(String.format("Finished in %s, speedup: %.2fx",
						stopwatch.getTime(2),
						(float)baseNs*numThreads/stopwatch.getTimeNs()
					));
				}
				
				// cleanup
				try {
					tasks.stopAndWait(10000);
				} catch (InterruptedException ex) {
					throw new Error(ex);
				}
			}
		}
	}
}
