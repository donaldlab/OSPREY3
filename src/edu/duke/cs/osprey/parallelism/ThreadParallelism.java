package edu.duke.cs.osprey.parallelism;

public class ThreadParallelism {

	private static int NUM_THREADS = Math.min(Runtime.getRuntime().availableProcessors(), 8);

	public static void setNumThreads( int threads ) { 
		NUM_THREADS = threads;

		if(NUM_THREADS < 1) NUM_THREADS = 1;

		else if(NUM_THREADS >= Runtime.getRuntime().availableProcessors()) 
			NUM_THREADS = Runtime.getRuntime().availableProcessors();
		
		setNumThreads();
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
