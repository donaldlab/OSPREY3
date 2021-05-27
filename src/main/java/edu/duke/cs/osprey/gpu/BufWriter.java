package edu.duke.cs.osprey.gpu;

import jdk.incubator.foreign.MemoryAddress;
import jdk.incubator.foreign.MemoryHandles;
import jdk.incubator.foreign.MemorySegment;

import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;


/**
 * A tool to easily write primitive types to a MemorySegment,
 * without the limitations of NIO ByteBuffers.
 */
public class BufWriter {

	private static final VarHandle hint = MemoryHandles.varHandle(int.class, ByteOrder.nativeOrder());
	private static final VarHandle hlong = MemoryHandles.varHandle(long.class, ByteOrder.nativeOrder());
	private static final VarHandle hfloat = MemoryHandles.varHandle(float.class, ByteOrder.nativeOrder());
	private static final VarHandle hdouble = MemoryHandles.varHandle(double.class, ByteOrder.nativeOrder());

	private static final long sizeInt = 4;
	private static final long sizeLong = 8;
	private static final long sizeFloat = 4;
	private static final long sizeDouble = 8;

	private final MemoryAddress base;

	public long pos = 0;

	public BufWriter(MemorySegment mem) {
		this.base = mem.address();
	}

	public void int32(long pos, int value) {
		hint.set(base.addOffset(pos), value);
	}
	public void int32(int value) {
		int32(pos, value);
		pos += sizeInt;
	}
	public long int32skip() {
		return skip(sizeInt);
	}

	public void uint32(long pos, long value) {
		hint.set(base.addOffset(pos), value & 0xffffffffL);
	}
	public void uint32(long value) {
		uint32(pos, value);
		pos += sizeInt;
	}
	public void uint32(long pos, int value) {
		hint.set(base.addOffset(pos), value);
	}
	public void uint32(int value) {
		uint32(pos, value);
		pos += sizeInt;
	}
	public long uint32skip() {
		return skip(sizeInt);
	}

	public void int64(long pos, long value) {
		hlong.set(base.addOffset(pos), value);
	}
	public void int64(long value) {
		int64(pos, value);
		pos += sizeLong;
	}
	public long int64skip() {
		return skip(sizeLong);
	}

	// sadly, Java can't represent uint64

	public void float32(long pos, float value) {
		hfloat.set(base.addOffset(pos), value);
	}
	public void float32(float value) {
		float32(pos, value);
		pos += sizeFloat;
	}
	public long float32skip() {
		return skip(sizeFloat);
	}

	public void float64(long pos, double value) {
		hdouble.set(base.addOffset(pos), value);
	}
	public void float64(double value) {
		float64(pos, value);
		pos += sizeDouble;
	}
	public long float64skip() {
		return skip(sizeDouble);
	}

	/**
	 * Skips the desired number of bytes,
	 * returns the position before skipping.
	 */
	public long skip(long bytes) {
		long pos = this.pos;
		this.pos += bytes;
		return pos;
	}

	public boolean isAligned(long alignment) {
		return pos % alignment == 0;
	}

	public void skipToAlignment(long alignment) {
		pos = padToAlignment(pos, alignment);
	}

	public static long padToAlignment(long size, long alignment) {
		// alignment won't be very big, so this is good enough for now
		while (size % alignment != 0) {
			size++;
		}
		return size;
	}

	public MemoryAddress address() {
		return base.addOffset(pos);
	}

	/**
	 * Sets the address of the header to the current position in the buffer,
	 * and increments the current position by the size of the header.
	 */
	public MemoryAddress place(Structs.Struct struct) {
		MemoryAddress addr = address();
		skip(struct.bytes());
		return addr;
	}

	public MemoryAddress place(Structs.Array array, long size) {
		MemoryAddress addr = address();
		skip(array.bytes(size));
		return addr;
	}

	public MemoryAddress place(Structs.Array array, long size, long align) {
		MemoryAddress addr = address();
		skip(array.bytes(size));
		skipToAlignment(align);
		return addr;
	}
}
