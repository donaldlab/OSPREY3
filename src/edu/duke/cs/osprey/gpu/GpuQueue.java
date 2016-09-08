package edu.duke.cs.osprey.gpu;

import com.jogamp.opencl.CLCommandQueue;

public class GpuQueue {
	
	private Gpu gpu;
	private CLCommandQueue queue;
	
	public GpuQueue(Gpu gpu, CLCommandQueue queue) {
		this.gpu = gpu;
		this.queue = queue;
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public CLCommandQueue getCLQueue() {
		return queue;
	}
	
	public void cleanup() {
		queue.release();
	}

	public boolean isProfilingEnabled() {
		return queue.isProfilingEnabled();
	}
}
