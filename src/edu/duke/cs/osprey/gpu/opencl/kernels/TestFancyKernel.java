package edu.duke.cs.osprey.gpu.opencl.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLMemory;

import edu.duke.cs.osprey.gpu.opencl.GpuQueue;
import edu.duke.cs.osprey.gpu.opencl.Kernel;

public class TestFancyKernel extends Kernel {
	
	private CLBuffer<DoubleBuffer> bufA;
	private CLBuffer<DoubleBuffer> bufB;
	private CLBuffer<DoubleBuffer> bufOut;
	
	private int workSize;
	private int groupSize;
		
	public TestFancyKernel(GpuQueue queue)
	throws IOException {
		super(queue, "test.cl", "fancy");
	}
	
	public DoubleBuffer getA() {
		return bufA.getBuffer();
	}
	
	public DoubleBuffer getB() {
		return bufB.getBuffer();
	}
	
	public DoubleBuffer getOut() {
		return bufOut.getBuffer();
	}
	
	public void setArgs(int workSize) {
		groupSize = getMaxGroupSize();
		this.workSize = roundUpWorkSize(workSize, groupSize);
		cleanup();
		CLContext context = getQueue().getCLQueue().getContext();
		bufA = context.createDoubleBuffer(workSize, CLMemory.Mem.READ_ONLY);
		bufB = context.createDoubleBuffer(workSize, CLMemory.Mem.READ_ONLY);
		bufOut = context.createDoubleBuffer(workSize, CLMemory.Mem.WRITE_ONLY);
	}
	
	public void runAsync() {
		getCLKernel()
			.putArg(bufA)
			.putArg(bufB)
			.putArg(bufOut)
			.rewind();
		runAsync(workSize, groupSize);
	}

	public void uploadSync() {
		uploadBufferAsync(bufA);
		uploadBufferAsync(bufB);
		waitForGpu();
	}

	public void downloadSync() {
		downloadBufferSync(bufOut);
	}
	
	public void cleanup() {
		if (bufA != null) {
			bufA.release();
		}
		if (bufB != null) {
			bufB.release();
		}
		if (bufOut != null) {
			bufOut.release();
		}
	}
}
