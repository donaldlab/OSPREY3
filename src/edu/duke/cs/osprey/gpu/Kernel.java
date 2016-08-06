package edu.duke.cs.osprey.gpu;

import java.io.IOException;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLKernel;

public abstract class Kernel<T extends BoundKernel<T>> {
	
	private CLKernel kernel;
	
	public Kernel(String filename, String name)
	throws IOException {
		kernel = Gpus.get().getProgram(filename).createCLKernel(name);
	}
	
	public CLKernel getCLKernel() {
		return kernel;
	}
	
	public abstract T bind(CLCommandQueue queue);
	
	public void cleanup() {
		kernel.release();
	}
}
