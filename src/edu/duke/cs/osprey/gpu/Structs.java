package edu.duke.cs.osprey.gpu;

import jdk.incubator.foreign.MemoryAddress;
import jdk.incubator.foreign.MemoryHandles;

import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;


/**
 * A crude way to represent and access C structs and arrays
 * using MemoryAddresses and MemoryHandles from the Foreign-Memory Access API.
 */
public class Structs {

	public static abstract class Struct {

		private Field[] fields = null;

		public <T extends Struct> T init(String ... fieldNames) {

			// get the fields
			Class<?> c = getClass();
			fields = new Field[fieldNames.length];
			for (int i=0; i<fieldNames.length; i++) {
				String fieldName = fieldNames[i];
				try {
					var declaredField = c.getDeclaredField(fieldName);
					declaredField.setAccessible(true);
					fields[i] = (Field)declaredField.get(this);
					fields[i].name = fieldName;
				} catch (NoSuchFieldException ex) {
					throw new Error("Can't initialize field: " + fieldName, ex);
				} catch (IllegalAccessException ex) {
					throw new Error("Can't read field: " + fieldName, ex);
				}
			}

			// find any missing fields
			Set<String> missingFields = Arrays.stream(c.getDeclaredFields())
				.map(f -> f.getName())
				.collect(Collectors.toSet());
			for (String name : fieldNames) {
				missingFields.remove(name);
			}
			missingFields.remove("this$0");
			if (!missingFields.isEmpty()) {
				throw new IllegalArgumentException("no order given for fields: " + missingFields);
			}

			@SuppressWarnings("unchecked")
			T struct = (T)this;
			return struct;
		}

		public void setAddress(MemoryAddress addr) {
			for (Field field : fields) {
				field.addr = addr;
				addr = addr.addOffset(field.bytes);
			}
		}

		/**
		 * Calculates the static size of the struct.
		 */
		public long bytes() {
			return Arrays.stream(fields)
				.mapToLong(Field::bytes)
				.sum();
		}

		/**
		 * Calculates the size of the struct assuming
		 * the last field is a dynamically-sized array.
		 */
		public long bytes(long arraySize) {
			long bytes = 0;
			for (int i=0; i<fields.length - 1; i++) {
				bytes += fields[i].bytes();
			}
			bytes += ((Array)fields[fields.length - 1]).bytes(arraySize);
			return bytes;
		}
	}

	public static abstract class Field {

		public static final long UnknownSize = -1;

		private final long bytes;

		protected String name;
		protected MemoryAddress addr;

		public Field(long bytes) {
			this.bytes = bytes;
		}

		public long bytes() {
			if (bytes == UnknownSize) {
				throw new IllegalStateException("field is dynamically sized");
			}
			return bytes;
		}
	}

	public static abstract class Array extends Field {

		public long staticSize;
		public final long itemBytes;

		public Array(long staticSize, long itemBytes) {
			super(staticSize == Field.UnknownSize ? Field.UnknownSize : staticSize*itemBytes);
			this.staticSize = staticSize;
			this.itemBytes = itemBytes;
		}

		public long bytes(long size) {
			if (size == Field.UnknownSize) {
				throw new IllegalArgumentException("size must be known");
			}
			return size*itemBytes;
		}

		public MemoryAddress offset(long i) {
			return addr.addOffset(i*itemBytes);
		}

		public void setAddress(MemoryAddress addr) {
			this.addr = addr;
		}
	}

	public static class Pad extends Field {

		public Pad(long bytes) {
			super(bytes);
		}
	}
	public static Pad pad(long bytes) {
		return new Pad(bytes);
	}

	public static class Int32 extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(int.class, ByteOrder.nativeOrder());

		public Int32() {
			super(4);
		}

		public int get() {
			return (int)handle.get(addr);
		}

		public void set(int value) {
			handle.set(addr, value);
		}
	}
	public static Int32 int32() {
		return new Int32();
	}

	public static class Uint32 extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(int.class, ByteOrder.nativeOrder());

		public Uint32() {
			super(4);
		}

		public long get() {
			return (int)handle.get(addr) & 0x00000000ffffffffL;
		}

		public void set(int value) {
			handle.set(addr, value);
		}

		public void set(long value) {
			handle.set(addr, (int)(value & 0x00000000ffffffffL));
		}
	}
	public static Uint32 uint32() {
		return new Uint32();
	}

	public static class Int64 extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(long.class, ByteOrder.nativeOrder());

		public Int64() {
			super(8);
		}

		public long get() {
			return (long)handle.get(addr);
		}

		public void set(long value) {
			handle.set(addr, value);
		}

		public static class Array extends Structs.Array {

			public Array(long size) {
				super(size, 8);
			}

			public long get(long i) {
				return (long)handle.get(offset(i));
			}

			public void set(long i, long value) {
				handle.set(offset(i), value);
			}
		}
	}
	public static Int64 int64() {
		return new Int64();
	}
	public static Int64.Array int64array(long size) {
		return new Int64.Array(size);
	}

	// sadly, Java can't represent a uint64


	public static class Float32 extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(float.class, ByteOrder.nativeOrder());

		public Float32() {
			super(4);
		}

		public float get() {
			return (float)handle.get(addr);
		}

		public void set(float value) {
			handle.set(addr, value);
		}

		public static class Array extends Structs.Array {

			public Array(long size) {
				super(size, 4);
			}

			public float get(long i) {
				return (float)handle.get(offset(i));
			}

			public void set(long i, float value) {
				handle.set(offset(i), value);
			}
		}
	}
	public static Float32 float32() {
		return new Float32();
	}
	public static Float32.Array float32array(long size) {
		return new Float32.Array(size);
	}

	public static class Float64 extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(double.class, ByteOrder.nativeOrder());

		public Float64() {
			super(8);
		}

		public double get() {
			return (double)handle.get(addr);
		}

		public void set(double value) {
			handle.set(addr, value);
		}

		public static class Array extends Structs.Array {

			public Array(long size) {
				super(size, 8);
			}

			public double get(long i) {
				return (double)handle.get(offset(i));
			}

			public void set(long i, double value) {
				handle.set(offset(i), value);
			}
		}
	}
	public static Float64 float64() {
		return new Float64();
	}
	public static Float64.Array float64array(long size) {
		return new Float64.Array(size);
	}

	// TODO: add more types
}
