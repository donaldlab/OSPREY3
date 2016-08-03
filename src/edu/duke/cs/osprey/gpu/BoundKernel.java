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
	private Gpu gpu;
	private CLCommandQueue queue;
	private CLEventList events;
	
	public BoundKernel(Kernel<T> kernel, Gpu gpu, boolean useProfiling) {
		this.kernel = kernel;
		this.gpu = gpu;
		if (useProfiling) {
			this.queue = gpu.getDevice().createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);
		} else {
			this.queue = gpu.getDevice().createCommandQueue();
		}
		this.events = null;
		
		if (this.queue.isOutOfOrderModeEnabled()) {
			throw new Error("GPU command queue should be strictly ordered... this is a bug");
		}
	}
	
	public Kernel<T> getKernel() {
		return kernel;
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public void waitForGpu() {
		queue.finish();
	}
	
	public void cleanup() {
		kernel.cleanup();
		queue.release();
	}
	
	protected int getMaxGroupSize() {
		return gpu.getDevice().getMaxWorkGroupSize();
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
		return gpu.getDevice().getContext().createBuffer(buf, type);
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
	
	public void initProfilingEvents() {
		
		// just in case...
		if (!queue.isProfilingEnabled()) {
			throw new IllegalStateException("profiling not enabled for this device. set Gpus.useProfiling = true before calling Gpus.get()");
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
