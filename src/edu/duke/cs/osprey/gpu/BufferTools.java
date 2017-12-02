package edu.duke.cs.osprey.gpu;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.nio.IntBuffer;
import java.nio.LongBuffer;

public class BufferTools {

	public static class OutOfDirectMemoryError extends VirtualMachineError {

		public OutOfDirectMemoryError(int bytes, OutOfMemoryError err) {
			super("can't allocate " + bytes + " bytes of direct memory", err);
		}
	}
	
	public static enum Type {
		
		Normal {
			
			@Override
			public ByteBuffer makeByte(int size) {
				return ByteBuffer.allocate(size);
			}
			
			@Override
			public DoubleBuffer makeDouble(int size) {
				return DoubleBuffer.allocate(size);
			}
			
			@Override
			public IntBuffer makeInt(int size) {
				return IntBuffer.allocate(size);
			}
			
			@Override
			public LongBuffer makeLong(int size) {
				return LongBuffer.allocate(size);
			}
		},
		Direct {
			
			@Override
			public ByteBuffer makeByte(int size) {
				try {
					return ByteBuffer.allocateDirect(size).order(ByteOrder.nativeOrder());
				} catch (OutOfMemoryError ex) {
					throw new OutOfDirectMemoryError(size, ex);
				}
			}
			
			@Override
			public DoubleBuffer makeDouble(int size) {
				return makeByte(size*Double.BYTES).asDoubleBuffer();
			}
			
			@Override
			public IntBuffer makeInt(int size) {
				return makeByte(size*Integer.BYTES).asIntBuffer();
			}
			
			@Override
			public LongBuffer makeLong(int size) {
				return makeByte(size*Long.BYTES).asLongBuffer();
			}
		};
		
		public abstract ByteBuffer makeByte(int size);
		public abstract DoubleBuffer makeDouble(int size);
		public abstract IntBuffer makeInt(int size);
		public abstract LongBuffer makeLong(int size);
	}
}
