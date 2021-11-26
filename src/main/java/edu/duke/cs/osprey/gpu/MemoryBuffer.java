package edu.duke.cs.osprey.gpu;

import jdk.incubator.foreign.MemoryHandles;
import jdk.incubator.foreign.MemorySegment;
import jdk.incubator.foreign.ResourceScope;

import java.lang.invoke.VarHandle;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;


/**
 * A tool to easily write primitive types to a MemorySegment,
 * without the limitations of NIO ByteBuffers and without exposing the incubating jdk.incubator.foreign module
 * to API clients.
 */
public class MemoryBuffer implements AutoCloseable {

	private static final VarHandle hint = MemoryHandles.varHandle(int.class, ByteOrder.nativeOrder());
	private static final VarHandle hlong = MemoryHandles.varHandle(long.class, ByteOrder.nativeOrder());
	private static final VarHandle hfloat = MemoryHandles.varHandle(float.class, ByteOrder.nativeOrder());
	private static final VarHandle hdouble = MemoryHandles.varHandle(double.class, ByteOrder.nativeOrder());
	private static final VarHandle hbyte = MemoryHandles.varHandle(byte.class, ByteOrder.nativeOrder());

	private static final long sizeInt = Integer.BYTES;
	private static final long sizeLong = Long.BYTES;
	private static final long sizeFloat = Float.BYTES;
	private static final long sizeDouble = Double.BYTES;

	private final MemorySegment base;
	private final ResourceScope sharedScope = ResourceScope.newSharedScope();
	private long pos = 0;

	private MemoryBuffer(MemorySegment mem) {
		this.base = mem;
	}

	public MemoryBuffer(long size) {
		this.base = MemorySegment.allocateNative(size, sharedScope);
	}

	public static MemoryBuffer ofByteBuffer(ByteBuffer buf) {
		return new MemoryBuffer(MemorySegment.ofByteBuffer(buf));
	}

	public MemoryBuffer sliceFrom(long offset) {
		return new MemoryBuffer(this.base.asSlice(offset));
	}

	public ByteBuffer asByteBuffer() {
		return base.asByteBuffer();
	}

	public long getPos() {
		return pos;
	}

	public void int32(long pos, int value) {
		hint.set(base, pos, value);
	}
	public void int32(int value) {
		int32(pos, value);
		pos += sizeInt;
	}
	public long int32skip() {
		return skip(sizeInt);
	}

	public void uint32(long pos, long value) {
		hint.set(base, pos, value & 0xffffffffL);
	}
	public void uint32(long value) {
		uint32(pos, value);
		pos += sizeInt;
	}
	public void uint32(long pos, int value) {
		hint.set(base, pos, value);
	}
	public void uint32(int value) {
		uint32(pos, value);
		pos += sizeInt;
	}
	public long uint32skip() {
		return skip(sizeInt);
	}

	public void int64(long pos, long value) {
		hlong.set(base, pos, value);
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
		hfloat.set(base, pos, value);
	}
	public void float32(float value) {
		float32(pos, value);
		pos += sizeFloat;
	}
	public long float32skip() {
		return skip(sizeFloat);
	}

	public void float64(long pos, double value) {
		hdouble.set(base, pos, value);
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

	/**
     * Returns a view of the MemoryBuffer at its current position and increments the position by the number of
	 * bytes the struct takes.
	 */
	public MemoryBuffer place(Structs.Struct struct) {
		var slice = base.asSlice(pos);
		skip(struct.bytes());
		return new MemoryBuffer(slice);
	}

	public MemoryBuffer place(Structs.Array array, long size) {
		var slice = base.asSlice(pos);
		skip(array.bytes(size));
		return new MemoryBuffer(slice);
	}

	public MemoryBuffer place(Structs.Array array, long size, long align) {
		var slice = base.asSlice(pos);
		skip(array.bytes(size));
		skipToAlignment(align);
		return new MemoryBuffer(slice);
	}

	private <T> void put(Structs.Field f, T value, VarHandle handle) {
		handle.set(base, f.offset, value);
	}

	private <T> void put(Structs.Array f, long idx, T value, VarHandle handle) {
		handle.set(base, idx * f.itemBytes, value);
	}

	@SuppressWarnings("unchecked")
	private <T> T get(Structs.Field f, VarHandle handle) {
		return (T) handle.get(base, f.offset);
	}

	@SuppressWarnings("unchecked")
	private <T> T get(Structs.Array a, long idx, VarHandle handle) {
		return (T) handle.get(base, idx * a.itemBytes);
	}

	public void putInt(Structs.Field f, int value) {
		put(f, value, hint);
	}

	public void putInt(Structs.Array f, long i, int value) {
		put(f, i, value, hint);
	}

	public int getInt(Structs.Field f) {
		return get(f, hint);
	}

	public int getInt(Structs.Array f, long i) {
		return get(f, i, hint);
	}

	public void putLong(Structs.Field f, long value) {
		put(f, value, hlong);
	}

	public void putLong(Structs.Array f, long i, long value) {
		put(f, i, value, hlong);
	}

	public long getLong(long offset) {
		return (long) hlong.get(base, offset);
	}

	public long getLong(Structs.Field f) {
		return get(f, hlong);
	}

	public long getLong(Structs.Array f, long i) {
		return get(f, i, hlong);
	}

	public void putFloat(Structs.Field f, float value) {
		put(f, value, hfloat);
	}

	public void putFloat(Structs.Array f, long i, float value) {
		put(f, i, value, hfloat);
	}

	public float getFloat(Structs.Field f) {
		return get(f, hfloat);
	}

	public float getFloat(Structs.Array f, long i) {
		return get(f, i, hfloat);
	}

	public void putDouble(Structs.Field f, double value) {
		put(f, value, hdouble);
	}

	public void putDouble(Structs.Array f, long i, double value) {
		put(f, i, value, hdouble);
	}

	public double getDouble(Structs.Field f) {
		return get(f, hdouble);
	}

	public double getDouble(Structs.Array f, long i) {
		return get(f, i, hdouble);
	}

	public void putBoolean(Structs.Field f, boolean value) {
		put(f, value ? (byte)1 : (byte)0, hbyte);
	}

	public void putBoolean(Structs.Array f, long i, boolean value) {
		put(f, i, value ? (byte)1 : (byte)0, hbyte);
	}

	public boolean getBoolean(Structs.Field f) {
		byte b = get(f, hbyte);
		return b != 0;
	}

	public boolean getBoolean(Structs.Array f, long i) {
		byte b = get(f, i, hbyte);
		return b != 0;
	}

	public void putChar(Structs.Field f, char value) {
		put(f, value, hbyte);
	}

	public void putChar(Structs.Array f, long i,  char value) {
		put(f, i, value, hbyte);
	}

	public char getChar(Structs.Field f) {
		return get(f, hbyte);
	}

	public char getChar(Structs.Array f, long i) {
		return get(f, i, hbyte);
	}

	public void putByte(Structs.Field f, byte value) {
		put(f, value, hbyte);
	}

	public void putByte(Structs.Array f, long i,  byte value) {
		put(f, i, value, hbyte);
	}

	public byte getByte(Structs.Field f) {
		return get(f, hbyte);
	}

	public byte getByte(Structs.Array f, long i) {
		return get(f, i, hbyte);
	}

	@Override
	public void close() throws IllegalStateException {
		sharedScope.close();
	}
}
