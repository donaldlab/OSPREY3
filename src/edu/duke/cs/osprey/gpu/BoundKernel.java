package edu.duke.cs.osprey.gpu;

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
	
	public void setWorkSize(int workSize) {
		this.workSize = roundUpWorkSize(workSize, gpu);
		initBuffers(workSize);
	}
	
	protected abstract void initBuffers(int workSize);
	public abstract void uploadAsync();
	public abstract void downloadSync();
	
	public void uploadSync() {
		uploadAsync();
		waitForGpu();
	}
	
	public void runAsync() {
		gpu.getQueue().put1DRangeKernel(kernel.getCLKernel(), 0, workSize, groupSize);
	}
	
	public void runSync() {
		runAsync();
		waitForGpu();
	}
	
	public void uploadRunDownloadSync() {
		uploadAsync();
		runAsync();
		waitForGpu();
		downloadSync();
	}
	
	public void waitForGpu() {
		gpu.getQueue().finish();
	}
	
	private int roundUpWorkSize(int workSize, Gpu gpu) {
		int groupSize = gpu.getDevice().getMaxWorkGroupSize();
		int r = workSize % groupSize;
		if (r == 0) {
			return workSize;
		} else {
			return (workSize + groupSize - 1)/groupSize*groupSize;
		}
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
	
	protected void uploadBufferAsync(CLBuffer<?> buf) {
		getGpu().getQueue().putWriteBuffer(buf, false);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		getGpu().getQueue().putReadBuffer(buf, true);
	}
}
