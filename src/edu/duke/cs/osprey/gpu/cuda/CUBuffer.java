package edu.duke.cs.osprey.gpu.cuda;

import java.nio.Buffer;

import com.jogamp.common.nio.Buffers;

import jcuda.Pointer;
import jcuda.driver.CUdeviceptr;

public class CUBuffer<T extends Buffer> {
	
	private Context context;
	private T buf;
	private long numBytes;
	private Pointer phBuf;
	private CUdeviceptr pdBuf;
	
	public CUBuffer(Context context, T buf) {
		
		this.context = context;
		this.buf = buf;
		
		numBytes = buf.capacity()*Buffers.sizeOfBufferElem(buf);
		
		// make the host pointer
		phBuf = Pointer.to(buf);
		
		// allocate device buffer
		pdBuf = context.malloc(numBytes);
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
		context.uploadSync(pdBuf, phBuf, numBytes);
	}
	
	public void uploadAsync() {
		buf.rewind();
		context.uploadAsync(pdBuf, phBuf, numBytes);
	}
	
	public T downloadSync() {
		buf.rewind();
		context.downloadSync(phBuf, pdBuf, numBytes);
		return buf;
	}
	
	public void cleanup() {
		context.free(pdBuf);
	}
}
