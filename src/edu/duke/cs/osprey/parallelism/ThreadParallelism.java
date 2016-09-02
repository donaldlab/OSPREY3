package edu.duke.cs.osprey.parallelism;

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
