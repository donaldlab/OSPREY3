package edu.duke.cs.osprey.gpu;

import java.io.IOException;
import java.nio.Buffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLMemory;

public abstract class Kernel {
	
	private GpuQueue queue;
	private CLKernel kernel;
	private CLEventList events;
	
	public Kernel(GpuQueue queue, String fileName, String kernelName)
	throws IOException {
		
		if (queue.getCLQueue().isOutOfOrderModeEnabled()) {
			throw new Error("GPU command queue should be strictly ordered... this is a bug");
		}
		
		this.queue = queue;
		this.kernel = queue.getGpu().getProgram(fileName).createCLKernel(kernelName); 
		this.events = null;
	}
	
	public GpuQueue getQueue() {
		return queue;
	}
	
	protected CLKernel getCLKernel() {
		return kernel;
	}
	
	public void waitForGpu() {
		queue.getCLQueue().finish();
	}
	
	public void cleanup() {
		kernel.release();
		if (events != null) {
			events.release();
		}
	}
	
	protected int getMaxGroupSize() {
		return queue.getCLQueue().getDevice().getMaxWorkGroupSize();
	}
	
	protected int roundUpWorkSize(int workSize, int groupSize) {
		return (workSize + groupSize - 1)/groupSize*groupSize;
	}
	
	protected <B extends Buffer> CLBuffer<B> wrapBuffer(B buf, CLMemory.Mem type) {
		return queue.getCLQueue().getContext().createBuffer(buf, type);
	}
	
	protected void uploadBufferAsync(CLBuffer<?> buf) {
		queue.getCLQueue().putWriteBuffer(buf, false, events);
	}
	
	protected void runAsync(int workSize, int groupSize) {
		queue.getCLQueue().put1DRangeKernel(kernel, 0, workSize, groupSize, events);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, true, events);
	}
	
	protected void downloadBufferAsync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, false);
	}
	
	public void initProfilingEvents(int eventListSize) {
		
		// just in case...
		if (!queue.isProfilingEnabled()) {
			throw new IllegalStateException("profiling not enabled for the bound command queue");
		}
		
		if (events != null) {
			events.release();
		}
		
		events = new CLEventList(eventListSize);
	}
	
	public CLEventList getProfilingEvents() {
		return events;
	}
	
	public void clearProfilingEvents() {
		if (events != null) {
			events.release();
			events = null;
		}
	}
}
