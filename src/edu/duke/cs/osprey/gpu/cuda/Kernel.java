package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;

import jcuda.Pointer;
import jcuda.driver.CUfunction;
import jcuda.driver.CUmodule;
import jcuda.driver.JCudaDriver;

public class Kernel {
	
	public class Function {
		
		public int numBlocks;
		public int blockThreads;
		public int sharedMemBytes;
		
		private CUfunction func;
		private Pointer pArgs;
		
		public Function(String name) {
			func = new CUfunction();
			JCudaDriver.cuModuleGetFunction(func, module, name);
			pArgs = null;
			numBlocks = 1;
			blockThreads = 1;
			sharedMemBytes = 0;
		}
		
		public void setArgs(Pointer val) {
			pArgs = val;
		}
		
		public void runAsync() {
			getContext().launchKernel(func, numBlocks, blockThreads, sharedMemBytes, pArgs, stream);
		}
	}
	
	private GpuStream stream;
	private CUmodule module;
	
	public Kernel(GpuStream stream, String filename)
	throws IOException {
		
		if (stream == null) {
			throw new IllegalArgumentException("stream can't be null");
		}
		
		this.stream = stream;
		
		module = getContext().getKernel(filename);
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
	
	protected static int divUp(int num, int denom) {
		return (num + denom - 1)/denom;
	}
}
