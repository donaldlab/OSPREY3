package edu.duke.cs.osprey.gpu;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;

public class BufferTools {
	
	public static enum Type {
		
		Normal {
			@Override
			public ByteBuffer make(int size) {
				return ByteBuffer.allocate(size);
			}
		},
		Direct {
			@Override
			public ByteBuffer make(int size) {
				return ByteBuffer.allocateDirect(size).order(ByteOrder.nativeOrder());
			}
		},
		CudaPinned {
			@Override
			public ByteBuffer make(int size) {
				// TODO
				throw new Error("not implemented yet");
			}
		};
		
		public abstract ByteBuffer make(int size);
	}
	
	public static ByteBuffer makeByte(int size, Type type) {
		return ByteBuffer.allocateDirect(size).order(ByteOrder.nativeOrder());
	}
	
	public static DoubleBuffer makeDouble(int size, Type type) {
		return makeByte(size*Double.BYTES, type).asDoubleBuffer();
	}
	
	public static IntBuffer makeInt(int size, Type type) {
		return makeByte(size*Integer.BYTES, type).asIntBuffer();
	}
	
	public static LongBuffer makeLong(int size, Type type) {
		return makeByte(size*Long.BYTES, type).asLongBuffer();
	}
}
