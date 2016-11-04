package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.GpuStream;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import jcuda.Pointer;

public class TestDPKernel extends Kernel {
	
	private CUBuffer<DoubleBuffer> a;
	private CUBuffer<DoubleBuffer> b;
	private CUBuffer<DoubleBuffer> out;
	
	private Pointer pArgs;
	
	public TestDPKernel(GpuStream stream, int numElements)
	throws IOException {
		super(stream, "testDP", "loop");
		
		a = stream.makeDoubleBuffer(numElements);
		b = stream.makeDoubleBuffer(numElements);
		out = stream.makeDoubleBuffer(numElements);
		
		pArgs = Pointer.to(
			Pointer.to(new int[] { numElements }),
			a.makeDevicePointer(),
			b.makeDevicePointer(),
			out.makeDevicePointer()
		);
	}
	
	public DoubleBuffer getA() {
		return a.getHostBuffer();
	}
	
	public DoubleBuffer getB() {
		return b.getHostBuffer();
	}
	
	public DoubleBuffer getOut() {
		return out.getHostBuffer();
	}
	
	public void uploadAsync() {
		a.uploadAsync();
		b.uploadAsync();
	}
	
	public void runAsync() {
		runAsync(1, 1, 0, pArgs);
	}
	
	public DoubleBuffer downloadSync() {
		out.downloadSync();
		return out.getHostBuffer();
	}
	
	public void cleanup() {
		a.cleanup();
		b.cleanup();
		out.cleanup();
	}
}
