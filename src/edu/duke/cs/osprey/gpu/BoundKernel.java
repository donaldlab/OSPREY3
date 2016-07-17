package edu.duke.cs.osprey.gpu;

import java.nio.Buffer;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLMemory;

public abstract class BoundKernel<T extends BoundKernel<T>> {
	
	private Kernel<T> kernel;
	private Gpu gpu;
	private int groupSize;
	private int workSize;
	
	public BoundKernel(Kernel<T> kernel, Gpu gpu) {
		this.kernel = kernel;
		this.gpu = gpu;
		groupSize = gpu.getDevice().getMaxWorkGroupSize();
		workSize = 0;
	}
	
	public Kernel<T> getKernel() {
		return kernel;
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	protected void setWorkSize(int val) {
		workSize = val;
	}
	
	public void runAsync() {
		if (workSize <= 0) {
			throw new IllegalStateException("invalid work size: " + workSize);
		}
		gpu.getQueue().put1DRangeKernel(kernel.getCLKernel(), 0, workSize, groupSize);
	}
	
	public void runSync() {
		runAsync();
		waitForGpu();
	}
	
	public void waitForGpu() {
		gpu.getQueue().finish();
	}
	
	public void cleanup() {
		kernel.cleanup();
	}
	
	protected int roundUpWorkSize(int workSize) {
		
		// get the number of threads per work group
		int groupSize = gpu.getDevice().getMaxWorkGroupSize();
		
		// round up to the nearest multiple of groupSize
		return (workSize + groupSize - 1)/groupSize*groupSize;
	}
	
	protected int getNumGroups(int workSize) {
		return workSize/gpu.getDevice().getMaxWorkGroupSize();
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
		getGpu().getQueue().putWriteBuffer(buf, false);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		getGpu().getQueue().putReadBuffer(buf, true);
	}
}
