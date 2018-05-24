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

package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

import jcuda.CudaException;
import jcuda.Pointer;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class Kernel {
	
	public static interface SharedMemCalculator {
		
		int calcBytes(int blockThreads);
		
		public static class None implements SharedMemCalculator {
			@Override
			public int calcBytes(int blockThreads) {
				return 0;
			}
		}
	}
		
	public class Function {
		
		private CUfunction func;
		private Pointer pArgs;
		
		public int numBlocks;
		public int blockThreads;
		public SharedMemCalculator sharedMemCalc;
		
		public Function(String name) {
			func = new CUfunction();
			JCudaDriver.cuModuleGetFunction(func, module, name);
			pArgs = null;
			numBlocks = 1;
			blockThreads = 1;
			sharedMemCalc = new SharedMemCalculator.None();
		}
		
		public void setArgs(Pointer val) {
			pArgs = val;
		}
		
		public void runAsync() {
			getContext().launchKernel(func, numBlocks, blockThreads, sharedMemCalc.calcBytes(blockThreads), pArgs, stream);
		}
		
		public int calcMaxBlockThreads() {
			
			// try the most threads, then step back until something works
			Gpu gpu = getContext().getGpu();
			for (int blockThreads = gpu.getMaxBlockThreads(); blockThreads > 0; blockThreads -= gpu.getWarpThreads()) {
				
				// can we launch at this size?
				if (canLaunch(blockThreads)) {
					return blockThreads;
				}
			}
		
			// we ran out of threads, this is really bad
			throw new Error("can't determine thread count for kernel launch, all thread counts failed");
		}
		
		public int getBestBlockThreads(AtomicInteger blockThreads) {
			return blockThreads.updateAndGet((int val) -> {
				if (val == -1) {
					val = calcMaxBlockThreads();
				}
				return val;
			});
		}
		
		private boolean canLaunch(int blockThreads) {
			try {
				
				// try to launch it (without changing internal state)
				int numBlocks = 1;
				getContext().launchKernel(func, numBlocks, blockThreads, sharedMemCalc.calcBytes(blockThreads), pArgs, stream);
				stream.waitForGpu();
				return true;
				
			} catch (CudaException ex) {
				if (ex.getMessage().equalsIgnoreCase("CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES")) {
					return false;
				} else {
					throw ex;
				}
			}
		}
	}
	
	private GpuStream stream;
	private CUmodule module;
	
	public Kernel(GpuStream stream, String filename) {
		
		if (stream == null) {
			throw new IllegalArgumentException("stream can't be null");
		}
		
		this.stream = stream;
		
		// make sure this thread is attached to this context
		stream.getContext().attachCurrentThread();
		
		try {
			module = getContext().getKernel(filename);
		} catch (IOException ex) {
			throw new Error("can't load Cuda kernel: " + filename, ex);
		}
	}
	
	public GpuStream getStream() {
		return stream;
	}
	
	public Context getContext() {
		return stream.getContext();
	}
	
	public Function makeFunction(String name) {
		return new Function(name);
	}
	
	public void waitForGpu() {
		getStream().waitForGpu();
	}
}
