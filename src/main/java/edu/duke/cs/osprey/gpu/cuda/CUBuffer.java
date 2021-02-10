/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gpu.cuda;

import java.nio.Buffer;

import com.jogamp.common.nio.Buffers;

import jcuda.CudaException;
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
		try {
			stream.getContext().pinBuffer(phBuf, numBytes);
		} catch (CudaException ex) {
			if (ex.getMessage().equals("CUDA_ERROR_HOST_MEMORY_ALREADY_REGISTERED")) {
				throw new Error("can't create new buffer... usually caused by failure to cleanup old buffers");
			} else {
				throw ex;
			}
		}
		
		// allocate device buffer
		pdBuf = stream.getContext().malloc(numBytes);
	}
	
	public T getHostBuffer() {
		return buf;
	}
	
	public Pointer getDevicePointer() {
		return Pointer.to(pdBuf);
	}
	
	public int size() {
		return buf.capacity();
	}
	
	public long getNumBytes() {
		return numBytes;
	}
	
	public void uploadAsync() {
		buf.rewind();
		stream.getContext().uploadAsync(pdBuf, phBuf, numBytes, stream);
	}
	
	public void downloadAsync() {
		buf.rewind();
		stream.getContext().downloadAsync(phBuf, pdBuf, numBytes, stream);
	}
	
	public T downloadSync() {
		downloadAsync();
		stream.waitForGpu();
		return buf;
	}
	
	public void expand(T buf) {
		
		// is the existing buffer big enough?
		int newNumBytes = buf.capacity()*Buffers.sizeOfBufferElem(buf);
		if (newNumBytes <= numBytes) {
			
			// yup
			return;
		}
		
		// nope, resize host side
		stream.getContext().unpinBuffer(phBuf);
		this.buf = buf;
		phBuf = Pointer.to(this.buf);
		stream.getContext().pinBuffer(phBuf, newNumBytes);
		
		// resize the device side
		stream.getContext().free(pdBuf);
		pdBuf = stream.getContext().malloc(newNumBytes);
		numBytes = newNumBytes;
	}
	
	public void cleanup() {
		try {
			stream.getContext().attachCurrentThread();
			stream.getContext().unpinBuffer(phBuf);
			stream.getContext().free(pdBuf);
			phBuf = null;
			pdBuf = null;
		} catch (Throwable t) {
			t.printStackTrace(System.err);
		}
	}
	
	@Override
	protected void finalize()
	throws Throwable {
		try {
			if (phBuf != null || pdBuf != null) {
				System.err.println("WARNING: " + getClass().getName() + " was garbage collected, but not cleaned up. Attempting cleanup now");
				cleanup();
			}
		} finally {
			super.finalize();
		}
	}
}
