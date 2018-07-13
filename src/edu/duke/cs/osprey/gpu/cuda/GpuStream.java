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
import java.nio.ByteBuffer;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;

import edu.duke.cs.osprey.gpu.BufferTools;
import jcuda.driver.CUstream;
import jcuda.driver.JCudaDriver;

public class GpuStream {
	
	public final BufferPool<ByteBuffer> byteBuffers;
	public final BufferPool<IntBuffer> intBuffers;
	public final BufferPool<LongBuffer> longBuffers;
	public final BufferPool<DoubleBuffer> doubleBuffers;
	
	private Context context;
	private CUstream stream;
	
	public GpuStream(Context context) {
		this.context = context;
		
		stream = new CUstream();
		JCudaDriver.cuStreamCreate(stream, 0);
		
		byteBuffers = new BufferPool<ByteBuffer>(
			(Integer size) -> makeByteBuffer(size),
			(CUBuffer<ByteBuffer> buf, int size) -> makeOrExpandByteBuffer(buf, size)
		);
		intBuffers = new BufferPool<IntBuffer>(
			(Integer size) -> makeIntBuffer(size),
			(CUBuffer<IntBuffer> buf, int size) -> makeOrExpandIntBuffer(buf, size)
		);
		longBuffers = new BufferPool<LongBuffer>(
			(Integer size) -> makeLongBuffer(size),
			(CUBuffer<LongBuffer> buf, int size) -> makeOrExpandLongBuffer(buf, size)
		);
		doubleBuffers = new BufferPool<DoubleBuffer>(
			(Integer size) -> makeDoubleBuffer(size),
			(CUBuffer<DoubleBuffer> buf, int size) -> makeOrExpandDoubleBuffer(buf, size)
		);
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
		return new CUBuffer<>(this, BufferTools.Type.Direct.makeByte(size));
	}
	
	public CUBuffer<IntBuffer> makeIntBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.Type.Direct.makeInt(size));
	}
	
	public CUBuffer<DoubleBuffer> makeDoubleBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.Type.Direct.makeDouble(size));
	}
	
	public CUBuffer<LongBuffer> makeLongBuffer(int size) {
		return new CUBuffer<>(this, BufferTools.Type.Direct.makeLong(size));
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
		if (buf.size() >= size) {
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
			byteBuffers.cleanup();
			intBuffers.cleanup();
			longBuffers.cleanup();
			doubleBuffers.cleanup();
			try {
				JCudaDriver.cuStreamDestroy(stream);
			} catch (Throwable t) {
				t.printStackTrace(System.err);
			}
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
