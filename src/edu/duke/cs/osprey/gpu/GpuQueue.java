package edu.duke.cs.osprey.gpu;

import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;

public class GpuQueue {
	
	private Gpu gpu;
	private CLCommandQueue queue;
	private CLContext separateContext;
	
	public GpuQueue(Gpu gpu) {
		this(gpu, false, false);
	}
	
	public GpuQueue(Gpu gpu, boolean useProfiling, boolean makeSeparateContext) {
		
		if (makeSeparateContext) {
			separateContext = CLContext.create(gpu.getDevice());
			this.gpu = new Gpu(separateContext.getDevices()[0]);
		} else {
			separateContext = null;
			this.gpu = gpu;
		}
		
		CLDevice device = this.gpu.getDevice();
		if (useProfiling) {
			queue = device.createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);
		} else {
			queue = device.createCommandQueue();
		}
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public CLCommandQueue getCLQueue() {
		return queue;
	}
	
	public void cleanup() {
		queue.release();
		
		if (separateContext != null) {
			separateContext.release();
			separateContext = null;
		}
	}

	public boolean isProfilingEnabled() {
		return queue.isProfilingEnabled();
	}
}
