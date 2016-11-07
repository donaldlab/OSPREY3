package edu.duke.cs.osprey.gpu.cuda;

import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;

import edu.duke.cs.osprey.gpu.BufferTools;
import jcuda.driver.CUstream;
import jcuda.driver.JCudaDriver;

public class GpuStream {
	
	private Context context;
	private CUstream stream;
	
	public GpuStream(Context context) {
		this.context = context;
		
		stream = new CUstream();
		JCudaDriver.cuStreamCreate(stream, 0);
	}
	
	public Context getContext() {
		return context;
	}

	public CUstream getStream() {
		return stream;
	}
	
	public <T extends Buffer> CUBuffer<T> makeBuffer(T buf) {
		return new CUBuffer<T>(this, buf);
	}
	
	public CUBuffer<ByteBuffer> makeByteBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.makeByte(size, BufferTools.Type.Direct));
	}
	
	public CUBuffer<IntBuffer> makeIntBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.makeInt(size, BufferTools.Type.Direct));
	}
	
	public CUBuffer<DoubleBuffer> makeDoubleBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.makeDouble(size, BufferTools.Type.Direct));
	}
	
	public CUBuffer<LongBuffer> makeLongBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.makeLong(size, BufferTools.Type.Direct));
	}
	
	public void waitForGpu() {
		JCudaDriver.cuStreamSynchronize(stream);
	}
	
	public void cleanup() {
		JCudaDriver.cuStreamDestroy(stream);
		stream = null;
	}
}
