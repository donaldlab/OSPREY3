package edu.duke.cs.osprey.gpu;

import java.io.IOException;

import com.jogamp.opencl.CLKernel;

public abstract class Kernel<T extends BoundKernel<T>> {
	
	private Gpu gpu;
	private CLKernel kernel;
	
	public Kernel(Gpu gpu, String filename, String name)
	throws IOException {
		this.gpu = gpu;
		kernel = gpu.getProgram(filename).createCLKernel(name);
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public CLKernel getCLKernel() {
		return kernel;
	}
	
	public abstract T bind(GpuQueue queue);
	
	public void cleanup() {
		kernel.release();
	}
}
