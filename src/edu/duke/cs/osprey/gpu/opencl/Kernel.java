package edu.duke.cs.osprey.gpu.opencl;

import java.io.IOException;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLException;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLPlatform;
import com.jogamp.opencl.llb.CL;
import com.jogamp.opencl.util.CLUtil;

public abstract class Kernel {
	
	private GpuQueue queue;
	private CLKernel kernel;
	private ProfilingEvents events;
	
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
		kernel = null;
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (kernel != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}
	
	protected int getMaxGroupSize() {
		return queue.getCLQueue().getDevice().getMaxWorkGroupSize();
	}
	
	protected int roundUpWorkSize(int workSize, int groupSize) {
		return (workSize + groupSize - 1)/groupSize*groupSize;
	}
	
	protected void uploadBufferAsync(CLBuffer<?> buf) {
		queue.getCLQueue().putWriteBuffer(buf, false, events != null ? events.getCLEvents() : null);
	}
	
	protected void uploadPartialBufferAsync(CLBuffer<?> buf) {
		
		// only write the part of the buffer below the limit
		buf.getBuffer().rewind();
		int size = buf.getElementSize()*buf.getBuffer().limit();
		
		// use the lower-level interface to upload partial buffers
		CL cl = CLPlatform.getLowLevelCLInterface();
		final int ret = cl.clEnqueueWriteBuffer(
			queue.getCLQueue().getID(),
			buf.getID(),
			CLUtil.clBoolean(false),
			0,
			size,
			buf.getBuffer(),
			0,
			null,
			events != null ? events.getCLIds() : null
		);

		if(ret != CL.CL_SUCCESS) {
			throw CLException.newException(ret, "can not enqueue write-buffer: " + buf);
		}
	}
	
	protected void runAsync(int workSize, int groupSize) {
		queue.getCLQueue().put1DRangeKernel(kernel, 0, workSize, groupSize, events != null ? events.getCLEvents() : null);
	}
	
	protected void downloadBufferSync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, true, events != null ? events.getCLEvents() : null);
	}
	
	protected void downloadBufferAsync(CLBuffer<?> buf) {
		queue.getCLQueue().putReadBuffer(buf, false, events != null ? events.getCLEvents() : null);
	}
	
	public ProfilingEvents getProfilingEvents() {
		return events;
	}
	public void setProfilingEvents(ProfilingEvents val) {
		
		if (val == events) {
			return;
		}
		
		// just in case...
		if (!queue.isProfilingEnabled()) {
			throw new IllegalStateException("profiling not enabled for the gpu queue");
		}
		
		events = val;
	}
}
