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

import edu.duke.cs.osprey.parallelism.TaskExecutor.Task;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.TimeFormatter;

@SuppressWarnings("unused")
public class BenchmarkTaskExecutor {
	
	private static abstract class TaskFactory {
		
		public int numRuns;
		
		protected TaskFactory(int numRuns) {
			this.numRuns = numRuns;
		}
		
		public abstract Task<Void> makeTask();
	}
	
	public static void main(String[] args) {
	
		// shows best speedups we can get on this hardware (ie, near-zero sync over head)
		//benchmarkBigCpu();
		
		// shows speedups including sync overhead
		benchmarkManySmallCpu();
	}
	
	private static void benchmarkBigCpu() {
		benchmark(new TaskFactory(1) {
			@Override
			public Task<Void> makeTask() {
				return () -> {
					// spin a lot
					for (int i=0; i<1e8; i++);
					return null;
				};
			}
		});
	}
	
	private static void benchmarkManySmallCpu() {
		benchmark(new TaskFactory(10000) {
			@Override
			public Task<Void> makeTask() {
				return () -> {
					// spin a little
					for (int i=0; i<1e4; i++);
					return null;
				};
			}
		});
	}
	
	private static void benchmark(TaskFactory factory) {
		
		int[] numThreadsList = { 1, 2, 4 };//, 8, 16, 32 };
		int numTrials = 6;
		
		long baseNs = Long.MAX_VALUE;
		long totalNs = 0;
		for (int i=0; i<numTrials; i++) {
			
			System.out.print(String.format("Benchmarking main thread...  "));
			
			TaskExecutor tasks = new TaskExecutor();
			Stopwatch stopwatch = benchmark(factory, tasks, 0);
			
			baseNs = Math.min(baseNs, stopwatch.getTimeNs());
			if (i >= numTrials/2) {
				totalNs += stopwatch.getTimeNs();
			}
			System.out.println("Finished in " + stopwatch.getTime(2));
		}
		System.out.println("\tLast half avg time " + TimeFormatter.format(totalNs*2/numTrials, 2));
		
		for (int numThreads : numThreadsList) {
			totalNs = 0;
			for (int i=0; i<numTrials; i++) {
		
				System.out.print(String.format("Benchmarking %2d threads...  ", numThreads));
				
				ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
				tasks.queueSize = factory.numRuns;
				tasks.start(numThreads);
				
				Stopwatch stopwatch = benchmark(factory, tasks, baseNs);
				
				System.out.println(String.format("Finished in %s, speedup: %.2fx",
					stopwatch.getTime(2),
					(float)baseNs*numThreads/stopwatch.getTimeNs()
				));
				
				if (i >= numTrials/2) {
					totalNs += stopwatch.getTimeNs();
				}
				
				// cleanup
				tasks.stopAndWait(10000);
			}
			long avgNs = totalNs*2/numTrials;
			System.out.println(String.format("\tLast half avg time %s, speedup: %.2fx",
				TimeFormatter.format(avgNs, 2),
				(float)baseNs*numThreads/avgNs
			));
		}
	}
	
	private static Stopwatch benchmark(TaskFactory factory, TaskExecutor tasks, long baseNs) {
		
		// do the warmup
		for (int j=0; j<tasks.getParallelism(); j++) {
			tasks.submit(factory.makeTask(), (task) -> {});
		}
		tasks.waitForFinish();
		
		// time the tasks
		int numTasks = tasks.getParallelism()*factory.numRuns;
		Stopwatch stopwatch = new Stopwatch().start();
		for (int j=0; j<numTasks; j++) {
			tasks.submit(factory.makeTask(), (task) -> {});
		}
		tasks.waitForFinish();
		stopwatch.stop();
		
		return stopwatch;
	}
}
