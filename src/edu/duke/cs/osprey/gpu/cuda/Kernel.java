package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;

import jcuda.Pointer;
import jcuda.driver.CUfunction;
import jcuda.driver.JCudaDriver;

public class Kernel {
	
	private Context context;
	private CUfunction func;
	
	public Kernel(Context context, String filename, String funcName)
	throws IOException {
		
		if (context == null) {
			throw new IllegalArgumentException("context can't be null");
		}
		
		this.context = context;
		
		func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, context.getKernel(filename), funcName);
	}
	
	public Context getContext() {
		return context;
	}
	
	protected void runAsync(int numBlocks, int blockThreads, int sharedMemBytes, Pointer pArgs) {
		context.launchKernel(func, numBlocks, blockThreads, sharedMemBytes, pArgs);
	}
	
	public void waitForGpu() {
		context.waitForGpu();
	}
	
	protected static int calcNumBlocks(int numThreads, int blockThreads) {
		return (numThreads + blockThreads - 1)/blockThreads;
	}
}
