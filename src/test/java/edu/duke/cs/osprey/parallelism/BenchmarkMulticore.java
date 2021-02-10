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
