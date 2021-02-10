package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.tools.Stopwatch;


public enum Timing {

	/**
	 * Uniform 1 s batches of computation.
	 * Good for benchmarking when precise timing on short timescales is necessary.
	 */
	Precise {
		@Override
		public long workMs(Stopwatch stopwatch) {
			return 1000; // 1 s work batches
		}
	},

	/**
	 * Within the first 5 minutes of computation, use quicker batches, but gradually grow to 1 min batches.
	 * Good for long-running computations when we want to minimize overhead.
	 */
	Efficient {
		@Override
		public long workMs(Stopwatch stopwatch) {
			return Math.max(500L, (long)(Math.min(stopwatch.getTimeM()/5.0, 1.0)*60_000));
		}
	};

	public abstract long workMs(Stopwatch stopwatch);
}
