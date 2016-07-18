package edu.duke.cs.osprey.gpu;

import java.nio.Buffer;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLMemory;

public abstract class BoundKernel<T extends BoundKernel<T>> {
	
	private Kernel<T> kernel;
	private Gpu gpu;
	
	public BoundKernel(Kernel<T> kernel, Gpu gpu) {
		this.kernel = kernel;
		this.gpu = gpu;
	}
	
	public Kernel<T> getKernel() {
		return kernel;
	}
	
	public Gpu getGpu() {
		return gpu;
	}
	
	public void waitForGpu() {
		gpu.getQueue().finish();
	}
	
	public void cleanup() {
		kernel.cleanup();
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
		gpu.getQueue().putWriteBuffer(buf, false);
	}
	
	protected void runAsync(int workSize, int groupSize) {
		gpu.getQueue().put1DRangeKernel(kernel.getCLKernel(), 0, workSize, groupSize);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		gpu.getQueue().putReadBuffer(buf, true);
	}
}
