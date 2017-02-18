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
			cfp.params.getInt("MinimizationThreads", 1),
			cfp.params.getInt("MinimizationGpus", 0),
			cfp.params.getInt("MinimizationStreamsPerGpu", 1)
		);
	}
	
	/** the number of CPU threads to use */
	public final int numThreads;
	
	/** the number of GPUs to use */
	public final int numGpus;
	
	/** the number of parallel tasks to send to each GPU */
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
	
	/** get the maximum number of tasks to be be executed in parallel */
	public int getParallelism() {
		return type.getParallelism(this);
	}
	
	public TaskExecutor makeTaskExecutor(boolean isStreaming) {
		if (getParallelism() > 1) {
			ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor();
			tasks.start(getParallelism(), isStreaming ? 0 : 1);
			return tasks;
		} else {
			return new TaskExecutor();
		}
	}
	
	public void cleanupTaskExecutor(TaskExecutor tasks) {
		if (tasks instanceof ThreadPoolTaskExecutor) {
			((ThreadPoolTaskExecutor)tasks).cleanup();
		}
	}
}
