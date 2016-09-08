package edu.duke.cs.osprey.gpu;

import java.nio.Buffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLMemory;

public abstract class BoundKernel<T extends BoundKernel<T>> {
	
	public static final int EventListSize = 128;
	
	private Kernel<T> kernel;
	private GpuQueue queue;
	private CLEventList events;
	
	public BoundKernel(Kernel<T> kernel, GpuQueue queue) {
		this.kernel = kernel;
		this.queue = queue;
		this.events = null;
		
		if (this.queue.getCLQueue().isOutOfOrderModeEnabled()) {
			throw new Error("GPU command queue should be strictly ordered... this is a bug");
		}
	}
	
	public Kernel<T> getKernel() {
		return kernel;
	}
	
	public GpuQueue getQueue() {
		return queue;
	}
	
	public void waitForGpu() {
		queue.getCLQueue().finish();
	}
	
	public void cleanup() {
		kernel.cleanup();
		queue.getCLQueue().release();
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
		queue.getCLQueue().put1DRangeKernel(kernel.getCLKernel(), 0, workSize, groupSize, events);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, true, events);
	}
	
	protected void downloadBufferAsync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, false);
	}
	
	public void initProfilingEvents() {
		
		// just in case...
		if (!queue.isProfilingEnabled()) {
			throw new IllegalStateException("profiling not enabled for the bound command queue");
		}
		
		events = new CLEventList(EventListSize);
	}
	
	public CLEventList getProfilingEvents() {
		return events;
	}
	
	public void clearProfilingEvents() {
		events = null;
	}
}
