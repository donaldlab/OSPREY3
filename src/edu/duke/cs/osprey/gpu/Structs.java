package edu.duke.cs.osprey.gpu;

import jdk.incubator.foreign.MemoryAddress;
import jdk.incubator.foreign.MemoryHandles;

import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.Set;
import java.util.function.Function;
import java.util.function.ToLongFunction;
import java.util.stream.Collectors;


/**
 * A crude way to represent and access C structs and arrays
 * using MemoryAddresses and MemoryHandles from the Foreign-Memory Access API.
 */
public class Structs {

	public static abstract class Struct {

		private Field[] fields = null;

		public <T extends Struct> T init(int bytes, String ... fieldNames) {

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

			// check the field sizes
			if (bytes != bytes()) {
				throw new IllegalArgumentException("struct size (" + bytes() + ") is not expected size (" + bytes + ")");
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
			return sum(fields, f -> f.bytes);
		}
	}

	public static <T> long sum(T[] things, ToLongFunction<? super T> converter) {
		return Arrays.stream(things)
			.mapToLong(converter)
			.sum();
	}

	public static abstract class Field {

		private final long bytes;

		protected String name;
		protected MemoryAddress addr;

		public Field(long bytes) {
			this.bytes = bytes;
		}
	}

	public static abstract class Array {

		public final long itemBytes;
		protected MemoryAddress addr;

		public Array(long itemBytes) {
			this.itemBytes = itemBytes;
		}

		public long bytes(long size) {
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

			public Array() {
				super(8);
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
	public static Int64.Array int64array() {
		return new Int64.Array();
	}


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

			public Array() {
				super(4);
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
	public static Float32.Array float32array() {
		return new Float32.Array();
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

			public Array() {
				super(8);
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
	public static Float64.Array float64array() {
		return new Float64.Array();
	}

	public static class Bool extends Field {

		private static final VarHandle handle = MemoryHandles.varHandle(byte.class, ByteOrder.nativeOrder());

		public Bool() {
			super(1);
		}

		public boolean get() {
			return (byte)handle.get(addr) != 0;
		}

		public void set(boolean value) {
			handle.set(addr, value ? (byte)1 : (byte)0);
		}

		// TODO: need bool array?
	}
	public static Bool bool() {
		return new Bool();
	}
}
