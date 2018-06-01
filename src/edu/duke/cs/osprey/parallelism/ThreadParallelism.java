/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
