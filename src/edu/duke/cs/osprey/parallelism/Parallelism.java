package edu.duke.cs.osprey.parallelism;

import edu.duke.cs.osprey.control.ConfigFileParser;

public class Parallelism {
	
	public static enum Type {
		
		Cpu {
			@Override
			public int getParallelism(Parallelism parallelism) {
				return parallelism.numThreads;
			}
		},
		Gpu {
			@Override
			public int getParallelism(Parallelism parallelism) {
				return parallelism.numGpus*parallelism.numStreamsPerGpu;
			}
		};
		
		public abstract int getParallelism(Parallelism parallelism);
	}
	
	public static Parallelism makeDefault() {
		return makeCpu(1);
	}
	
	public static Parallelism makeCpu(int numThreads) {
		return new Parallelism(numThreads, 0, 0);
	}
	
	public static Parallelism makeGpu(int numGpus, int numStreamsPerGpu) {
		return new Parallelism(0, numGpus, numStreamsPerGpu);
	}
	
	// TODO: this should eventually go into a CFP-only area
	// it can be moved when we start refactoring config stuff to prepare for Python-land
	public static Parallelism makeFromConfig(ConfigFileParser cfp) {
		return new Parallelism(
			cfp.getParams().getInt("MinimizationThreads", 1),
			cfp.getParams().getInt("MinimizationGpus", 0),
			cfp.getParams().getInt("MinimizationStreamsPerGpu", 1)
		);
	}
	
	public final int numThreads;
	public final int numGpus;
	public final int numStreamsPerGpu;
	public final Type type;
	
	public Parallelism(int numThreads, int numGpus, int numStreamsPerGpu) {
		this.numThreads = numThreads;
		this.numGpus = numGpus;
		this.numStreamsPerGpu = numStreamsPerGpu;
		
		// prefer gpus over threads
		if (numGpus > 0) {
			type = Type.Gpu;
		} else {
			type = Type.Cpu;
		}
		
		if (getParallelism() <= 0) {
			throw new IllegalArgumentException(String.format("parallelism should be at least 1: threads=%d, gpus=%d, streams/gpu=%d",
				numThreads, numGpus, numStreamsPerGpu
			));
		}
	}
	
	public int getParallelism() {
		return type.getParallelism(this);
	}
	
	public TaskExecutor makeTaskExecutor() {
		return makeTaskExecutor(null);
	}
	
	/**
	 * Makes a TaskExecutor to process tasks, possibly in parallel
	 * @param useQueue true to buffer tasks in a queue before processing (can be faster),
	 *                 false to only submit a task when a thread is ready (prevents extra tasks)
	 */
	public TaskExecutor makeTaskExecutor(Integer queueSize) {
		if (getParallelism() > 1) {
			ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
			if (queueSize != null) {
				tasks.queueSize = queueSize;
			}
			tasks.start(getParallelism());
			return tasks;
		} else {
			return new TaskExecutor();
		}
	}
}
