package edu.duke.cs.osprey.gpu.cuda;

import java.nio.Buffer;

import com.jogamp.common.nio.Buffers;

import jcuda.Pointer;
import jcuda.driver.CUdeviceptr;

public class CUBuffer<T extends Buffer> {
	
	private GpuStream stream;
	private T buf;
	private long numBytes;
	private Pointer phBuf;
	private CUdeviceptr pdBuf;
	
	public CUBuffer(GpuStream stream, T buf) {
		
		this.stream = stream;
		this.buf = buf;
		
		numBytes = buf.capacity()*Buffers.sizeOfBufferElem(buf);
		
		// make the host pointer
		phBuf = Pointer.to(buf);
		
		// allocate device buffer
		pdBuf = stream.getContext().malloc(numBytes);
	}
	
	public T getHostBuffer() {
		return buf;
	}
	
	public Pointer makeDevicePointer() {
		return Pointer.to(pdBuf);
	}
	
	public long getNumBytes() {
		return numBytes;
	}
	
	public void uploadSync() {
		stream.getContext().uploadSync(pdBuf, phBuf, numBytes);
	}
	
	public void uploadAsync() {
		buf.rewind();
		stream.getContext().uploadAsync(pdBuf, phBuf, numBytes, stream);
	}
	
	public T downloadSync() {
		buf.rewind();
		stream.getContext().downloadSync(phBuf, pdBuf, numBytes);
		return buf;
	}
	
	public void cleanup() {
		stream.getContext().free(pdBuf);
	}
}
