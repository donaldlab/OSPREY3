package edu.duke.cs.osprey.gpu.cuda.kernels;

import java.io.IOException;
import java.nio.DoubleBuffer;

import edu.duke.cs.osprey.gpu.BufferTools;
import edu.duke.cs.osprey.gpu.cuda.CUBuffer;
import edu.duke.cs.osprey.gpu.cuda.Context;
import edu.duke.cs.osprey.gpu.cuda.Kernel;
import jcuda.Pointer;

public class TestKernel extends Kernel {
	
	private CUBuffer<DoubleBuffer> a;
	private CUBuffer<DoubleBuffer> b;
	private CUBuffer<DoubleBuffer> out;
	
	private int numBlocks;
	private int blockThreads;
	
	private Pointer pArgs;
	
	public TestKernel(Context context, int numElements)
	throws IOException {
		super(context, "test.cubin", "add");
		
		a = new CUBuffer<>(getContext(), BufferTools.makeDouble(numElements, BufferTools.Type.Direct));
		b = new CUBuffer<>(getContext(), BufferTools.makeDouble(numElements, BufferTools.Type.Direct));
		out = new CUBuffer<>(getContext(), BufferTools.makeDouble(numElements, BufferTools.Type.Direct));
		
		pArgs = Pointer.to(
			Pointer.to(new int[] { numElements }),
			a.makeDevicePointer(),
			b.makeDevicePointer(),
			out.makeDevicePointer()
		);
		
		blockThreads = 256;
		numBlocks = calcNumBlocks(numElements, blockThreads);
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
		//a.uploadAsync();
		//b.uploadAsync();
		a.uploadSync();
		b.uploadSync();
	}
	
	public void runAsync() {
		runAsync(numBlocks, blockThreads, 0, pArgs);
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
