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

package edu.duke.cs.osprey;

import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.parallelism.TimingThread;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.ArrayList;
import java.util.List;


public class Benchmark {

	public final Stopwatch stopwatch;
	public final long numOps;
	public final double opsPerSecond;

	public Benchmark(int numRuns, Runnable task) {
		this(1, numRuns, task);
	}

	public Benchmark(int numWarmups, int numRuns, Runnable task) {
		this(1, numWarmups, numRuns, task);
	}

	public Benchmark(int numThreads, int numWarmups, int numRuns, Runnable task) {

		// make the threads
		List<TimingThread> threads = new ArrayList<>();
		for (int i=0; i<numThreads; i++) {
			threads.add(new TimingThread("Benchmark-" + i) {

				private Molecule mol;
				private EnergyFunction efunc;

				@Override
				public void warmup() {
					for (int i=0; i<numWarmups; i++) {
						task.run();
					}
				}

				@Override
				public void time() {
					for (int i=0; i<numRuns; i++) {
						task.run();
					}
				}
			});
		}

		// run threads and wait
		stopwatch = TimingThread.timeThreads(threads);

		// compute the summary statistics
		numOps = numThreads*numRuns;
		opsPerSecond = numOps/stopwatch.getTimeS();
	}

	@Override
	public String toString() {
		return String.format("finished %d ops in %s @ %.2f ops/s",
			numOps,
			stopwatch.getTime(2),
			opsPerSecond
		);
	}

	public String toString(Benchmark other) {
		return toString() + String.format(", speedup %.2fx",
			this.opsPerSecond/other.opsPerSecond
		);
	}
}
