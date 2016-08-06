package edu.duke.cs.osprey.gpu;

import java.nio.Buffer;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLEventList;
import com.jogamp.opencl.CLMemory;

public abstract class BoundKernel<T extends BoundKernel<T>> {
	
	public static final int EventListSize = 128;
	
	private Kernel<T> kernel;
	private CLCommandQueue queue;
	private CLEventList events;
	
	public BoundKernel(Kernel<T> kernel, CLCommandQueue queue) {
		this.kernel = kernel;
		this.queue = queue;
		this.events = null;
		
		if (this.queue.isOutOfOrderModeEnabled()) {
			throw new Error("GPU command queue should be strictly ordered... this is a bug");
		}
	}
	
	public Kernel<T> getKernel() {
		return kernel;
	}
	
	public CLCommandQueue getQueue() {
		return queue;
	}
	
	public void waitForGpu() {
		queue.finish();
	}
	
	public void cleanup() {
		kernel.cleanup();
		queue.release();
	}
	
	protected int getMaxGroupSize() {
		return queue.getDevice().getMaxWorkGroupSize();
	}
	
	protected int roundUpWorkSize(int workSize, int groupSize) {
		return (workSize + groupSize - 1)/groupSize*groupSize;
	}
	
	protected CLBuffer<DoubleBuffer> makeOrIncreaseBuffer(CLBuffer<DoubleBuffer> buf, int workSize) {
		return makeOrIncreaseBuffer(buf, workSize, CLMemory.Mem.READ_WRITE);
	}
	
	protected CLBuffer<DoubleBuffer> makeOrIncreaseBuffer(CLBuffer<DoubleBuffer> buf, int workSize, CLMemory.Mem type) {
		if (buf == null || buf.getCLCapacity() < workSize) {
			buf = Gpus.get().getContext().createDoubleBuffer(workSize, type);
		}
		return buf;
	}
	
	protected <B extends Buffer> CLBuffer<B> wrapBuffer(B buf, CLMemory.Mem type) {
		return queue.getContext().createBuffer(buf, type);
	}
	
	protected void uploadBufferAsync(CLBuffer<?> buf) {
		queue.putWriteBuffer(buf, false, events);
	}
	
	protected void runAsync(int workSize, int groupSize) {
		queue.put1DRangeKernel(kernel.getCLKernel(), 0, workSize, groupSize, events);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		queue.putReadBuffer(buf, true, events);
	}
	
	protected void downloadBufferAsync(CLBuffer<?> buf) {
		queue.putReadBuffer(buf, false);
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
