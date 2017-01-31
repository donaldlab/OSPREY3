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
	
	public CUBuffer<ByteBuffer> makeOrExpandByteBuffer(CUBuffer<ByteBuffer> buf, int size) {
		buf = useOrCleanupBuffer(buf, size);
		if (buf == null) {
			buf = makeByteBuffer(size);
		}
		return buf;
	}
	
	public CUBuffer<IntBuffer> makeOrExpandIntBuffer(CUBuffer<IntBuffer> buf, int size) {
		buf = useOrCleanupBuffer(buf, size);
		if (buf == null) {
			buf = makeIntBuffer(size);
		}
		return buf;
	}
	
	public CUBuffer<DoubleBuffer> makeOrExpandDoubleBuffer(CUBuffer<DoubleBuffer> buf, int size) {
		buf = useOrCleanupBuffer(buf, size);
		if (buf == null) {
			buf = makeDoubleBuffer(size);
		}
		return buf;
	}
	
	public CUBuffer<LongBuffer> makeOrExpandLongBuffer(CUBuffer<LongBuffer> buf, int size) {
		buf = useOrCleanupBuffer(buf, size);
		if (buf == null) {
			buf = makeLongBuffer(size);
		}
		return buf;
	}
	
	private <T extends Buffer> CUBuffer<T> useOrCleanupBuffer(CUBuffer<T> buf, int size) {
		
		if (buf == null) {
			return null;
		}
		
		// if the old buffer is big enough, use that
		if (buf.getHostBuffer().capacity() >= size) {
			return buf;
		}
		
		// otherwise, clean it up
		buf.cleanup();
		return null;
	}
	
	public <T extends Buffer> CUBuffer<T> makeOrExpandBuffer(CUBuffer<T> dest, T src) {
		if (dest == null) {
			return makeBuffer(src);
		}
		dest.expand(src);
		return dest;
	}
	
	public void waitForGpu() {
		JCudaDriver.cuStreamSynchronize(stream);
	}
	
	public void cleanup() {
		if (stream != null) {
			JCudaDriver.cuStreamDestroy(stream);
			stream = null;
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (stream != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}
}
