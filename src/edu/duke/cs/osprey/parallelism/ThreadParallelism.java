package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.ForkJoinPool;

public class ThreadParallelism {

	// the default number of threads should probably always be 1
	// until the user decides to commit more resources to the process
	// especially since this value gets sent to MultiTermEnergyFunction,
	// which actually performs worse with more threads on small/medium designs!
	private static int NUM_THREADS = 1; //Runtime.getRuntime().availableProcessors()/2;
	
	static {
		// make sure the ForkJoinPool uses our default, and not its own default
		setForkJoinPoolConstructorArg(NUM_THREADS);
	}

	public static void setNumThreads( int threads ) {
		
		// I think we can only set this once per process
		// so at least warn a dev if we try to do that more than once
		boolean itWorked = setNumThreadsIfPossible(threads);
		assert (itWorked) : String.format("tried to set fork join pool to %d threads, but it has %d instead", NUM_THREADS, ForkJoinPool.commonPool().getParallelism());
	}
	
	public static boolean setNumThreadsIfPossible(int numThreads) {
		
		// short circuit if nothing changed
		if (numThreads == NUM_THREADS) {
			return true;
		}
		
		// clamp numThreads to acceptable range
		if (numThreads < 1) {
			numThreads = 1;
		} else if (numThreads >= Runtime.getRuntime().availableProcessors()) {
			numThreads = Runtime.getRuntime().availableProcessors();
		}
		
		// set the property that tells the automatically-constructor ForkJoinPool what to do
		// hopefully, it hasn't been constructed yet
		setForkJoinPoolConstructorArg(numThreads);
		
		// did it work?
		boolean itWorked = ForkJoinPool.commonPool().getParallelism() == numThreads;
		
		if (itWorked) {
			NUM_THREADS = numThreads;
		}
		
		return itWorked;
	}

	public static int getNumThreads() {
		return NUM_THREADS;
	}
	
	public static void setDefaultNumThreads() {
		setNumThreads( NUM_THREADS );
	}
	
	private static void setForkJoinPoolConstructorArg(int numThreads) {
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(numThreads));
	}
}
