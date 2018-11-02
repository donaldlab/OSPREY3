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
