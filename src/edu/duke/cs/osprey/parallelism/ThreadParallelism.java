package edu.duke.cs.osprey.parallelism;

import java.util.concurrent.ForkJoinPool;

public class ThreadParallelism {

	// the default number of threads should probably always be 1
	// until the user decides to commit more resources to the process
	// especially since this value gets sent to MultiTermEnergyFunction,
	// which actually performs worse with more threads on small/medium designs!
	private static int NUM_THREADS = 1; //Runtime.getRuntime().availableProcessors()/2;

	public static void setNumThreads( int threads ) { 
		NUM_THREADS = threads;

		if(NUM_THREADS < 1) NUM_THREADS = 1;

		else if(NUM_THREADS >= Runtime.getRuntime().availableProcessors()) 
			NUM_THREADS = Runtime.getRuntime().availableProcessors();
		
		setNumThreads();
		
		// I think we can only set this once per process
		// so at least warn a dev if we try to do that more than once
		assert (ForkJoinPool.commonPool().getParallelism() == NUM_THREADS)
			: String.format("tried to set fork join pool to %d threads, but it has %d instead", NUM_THREADS, ForkJoinPool.commonPool().getParallelism());
	}

	public static int getNumThreads() {
		return NUM_THREADS;
	}
	
	public static void setDefaultNumThreads() {
		setNumThreads( NUM_THREADS );
	}
	
	private static void setNumThreads() {
		System.setProperty( "java.util.concurrent.ForkJoinPool.common.parallelism", String.valueOf(getNumThreads()) );
	}

}
