package edu.duke.cs.osprey.gpu.cuda;

import java.io.IOException;

import jcuda.Pointer;
import jcuda.driver.CUfunction;
import jcuda.driver.JCudaDriver;

public class Kernel {
	
	private GpuStream stream;
	private CUfunction func;
	
	public Kernel(GpuStream stream, String filename, String funcName)
	throws IOException {
		
		if (stream == null) {
			throw new IllegalArgumentException("stream can't be null");
		}
		
		this.stream = stream;
		
		func = new CUfunction();
		JCudaDriver.cuModuleGetFunction(func, getContext().getKernel(filename), funcName);
	}
	
	public GpuStream getStream() {
		return stream;
	}
	
	public Context getContext() {
		return stream.getContext();
	}
	
	protected void runAsync(int numBlocks, int blockThreads, int sharedMemBytes, Pointer pArgs) {
		getContext().launchKernel(func, numBlocks, blockThreads, sharedMemBytes, pArgs, stream);
	}
	
	public void waitForGpu() {
		getStream().waitForGpu();
	}
	
	protected static int divUp(int num, int denom) {
		return (num + denom - 1)/denom;
	}
}
