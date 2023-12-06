package edu.duke.cs.osprey.gpu;

import java.lang.foreign.MemorySegment;
import java.lang.foreign.Arena;
import java.lang.foreign.ValueLayout;
import java.nio.ByteBuffer;


/**
 * A tool to easily write primitive types to a MemorySegment,
 * without the limitations of NIO ByteBuffers and without exposing the preview java.lang.foreign module
 * to API clients.
 * <p>
 * TODO: preview JEP 424, 434, and 442 probably obviate the purpose of this class. When java.lang.foreign is generally
 * available, consider removing this class
 */
public class MemoryBuffer implements AutoCloseable {
    private static final ValueLayout.OfInt hint = ValueLayout.JAVA_INT;
    private static final ValueLayout.OfLong hlong = ValueLayout.JAVA_LONG;
    private static final ValueLayout.OfFloat hfloat = ValueLayout.JAVA_FLOAT;
    private static final ValueLayout.OfDouble hdouble = ValueLayout.JAVA_DOUBLE;
    private static final ValueLayout.OfByte hbyte = ValueLayout.JAVA_BYTE;

    private static final long sizeInt = Integer.BYTES;
    private static final long sizeLong = Long.BYTES;
    private static final long sizeFloat = Float.BYTES;
    private static final long sizeDouble = Double.BYTES;

    private final MemorySegment base;
    private final Arena sharedScope = Arena.ofShared();
    private long pos = 0;

    private MemoryBuffer(MemorySegment mem) {
        this.base = mem;
    }

    public MemoryBuffer(long size) {
        this.base = sharedScope.allocate(size, 8);
    }

    public static MemoryBuffer ofByteBuffer(ByteBuffer buf) {
        return new MemoryBuffer(MemorySegment.ofBuffer(buf));
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
        base.set(hint, pos, value);
    }

    public void int32(int value) {
        int32(pos, value);
        pos += sizeInt;
    }

    public long int32skip() {
        return skip(sizeInt);
    }

    public void uint32(long pos, long value) {
        base.set(hint, pos, (int) (value & 0xffffffffL));
    }

    public void uint32(long value) {
        uint32(pos, value);
        pos += sizeInt;
    }

    public void uint32(long pos, int value) {
        base.set(hint, pos, value);
    }

    public void uint32(int value) {
        uint32(pos, value);
        pos += sizeInt;
    }

    public long uint32skip() {
        return skip(sizeInt);
    }

    public void int64(long pos, long value) {
        base.set(hlong, pos, value);
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
        base.set(hfloat, pos, value);
    }

    public void float32(float value) {
        float32(pos, value);
        pos += sizeFloat;
    }

    public long float32skip() {
        return skip(sizeFloat);
    }

    public void float64(long pos, double value) {
        base.set(hdouble, pos, value);
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

	/*
	private <T> void put(Structs.Field f, T value, ValueLayout handle) {
		handle.set(base, f.offset, value);
	}

	private <T> void put(Structs.Array f, long idx, T value, ValueLayout handle) {
		handle.set(base, idx * f.itemBytes, value);
	}

	@SuppressWarnings("unchecked")
	private <T> T get(Structs.Field f, ValueLayout handle) {
		return (T) handle.get(base, f.offset);
	}

	@SuppressWarnings("unchecked")
	private <T> T get(Structs.Array a, long idx, ValueLayout handle) {
		return (T) handle.get(base, idx * a.itemBytes);
	}
	 */

    public void putInt(Structs.Field f, int value) {
        base.set(hint, f.offset, value);
    }

    public void putInt(Structs.Array f, long i, int value) {
        base.set(hint, i * f.itemBytes, value);
    }

    public int getInt(Structs.Field f) {
        return base.get(hint, f.offset);
    }

    public int getInt(Structs.Array f, long i) {
        return base.get(hint, i * f.itemBytes);
    }

    public void putLong(Structs.Field f, long value) {
        base.set(hlong, f.offset, value);
    }

    public void putLong(Structs.Array f, long i, long value) {
        base.set(hlong, i * f.itemBytes, value);
    }

    public long getLong(long offset) {
        return base.get(hlong, offset);
    }

    public long getLong(Structs.Field f) {
        return base.get(hlong, f.offset);
    }

    public long getLong(Structs.Array f, long i) {
        return base.get(hlong, i * f.itemBytes);
    }

    public void putFloat(Structs.Field f, float value) {
        base.set(hfloat, f.offset, value);
    }

    public void putFloat(Structs.Array f, long i, float value) {
        base.set(hfloat, i * f.itemBytes, value);
    }

    public float getFloat(Structs.Field f) {
        return base.get(hfloat, f.offset);
    }

    public float getFloat(Structs.Array f, long i) {
        return base.get(hfloat, i * f.itemBytes);
    }

    public void putDouble(Structs.Field f, double value) {
        base.set(hdouble, f.offset, value);
    }

    public void putDouble(Structs.Array f, long i, double value) {
        base.set(hdouble, i * f.itemBytes, value);
    }

    public double getDouble(Structs.Field f) {
        return base.get(hdouble, f.offset);
    }

    public double getDouble(Structs.Array f, long i) {
        return base.get(hdouble, i * f.itemBytes);
    }

    public void putBoolean(Structs.Field f, boolean value) {
        base.set(hbyte, f.offset, value ? (byte) 1 : (byte) 0);
    }

    public void putBoolean(Structs.Array f, long i, boolean value) {
        base.set(hbyte, i * f.itemBytes, value ? (byte) 1 : (byte) 0);
    }

    public boolean getBoolean(Structs.Field f) {
        return base.get(hbyte, f.offset) != 0;
    }

    public boolean getBoolean(Structs.Array f, long i) {
        return base.get(hbyte, i * f.itemBytes) != 0;
    }

    public void putChar(Structs.Field f, char value) {
        base.set(hbyte, f.offset, (byte) value);
    }

    public void putChar(Structs.Array f, long i, char value) {
        base.set(hbyte, i * f.itemBytes, (byte) value);
    }

    public char getChar(Structs.Field f) {
        return (char) base.get(hbyte, f.offset);
    }

    public char getChar(Structs.Array f, long i) {
        return (char) base.get(hbyte, i * f.itemBytes);
    }

    public void putByte(Structs.Field f, byte value) {
        base.set(hbyte, f.offset, value);
    }

    public void putByte(Structs.Array f, long i, byte value) {
        base.set(hbyte, i * f.itemBytes, value);
    }

    public byte getByte(Structs.Field f) {
        return base.get(hbyte, f.offset);
    }

    public byte getByte(Structs.Array f, long i) {
        return base.get(hbyte, i * f.itemBytes);
    }

    @Override
    public void close() throws IllegalStateException {
        sharedScope.close();
    }
}
