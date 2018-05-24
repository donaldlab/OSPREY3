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

import edu.duke.cs.osprey.control.ConfigFileParser;

/** Specified how Osprey should use available hardware. */
public class Parallelism {
	
	public static class Builder {
		
		/** The number of CPUs cores to use. */
		private int numCpus = 1;
		
		/** The number of GPUs to use. */
		private int numGpus = 0;
		
		/** The number of simultaneous tasks that should be given to each GPU */
		private int numStreamsPerGpu = 1;
		
		public Builder setNumCpus(int val) {
			numCpus = val;
			return this;
		}
		
		public Builder setNumGpus(int val) {
			numGpus = val;
			return this;
		}
		
		public Builder setNumStreamsPerGpu(int val) {
			numStreamsPerGpu = val;
			return this;
		}
		
		public Parallelism build() {
			return new Parallelism(numCpus, numGpus, numStreamsPerGpu);
		}
	}
	
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
	
	public static Parallelism makeCpu(int numThreads) {
		return new Parallelism(numThreads, 0, 0);
	}
	
	public static int getMaxNumCPUs() {
		return Runtime.getRuntime().availableProcessors();
	}
	
	public static Parallelism make(int numCpus, int numGpus, int numStreamsPerGpu) {
		return new Parallelism(numCpus, numGpus, numStreamsPerGpu);
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
	
	/** get the maximum number of tasks to be be executed in parallel */
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
