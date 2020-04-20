package edu.duke.cs.osprey.energy.compiled;

import com.sun.jna.*;

import static edu.duke.cs.osprey.gpu.Structs.*;
import static edu.duke.cs.osprey.tools.Log.log;

import java.lang.invoke.VarHandle;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.gpu.BufWriter;
import edu.duke.cs.osprey.gpu.Structs;
import jdk.incubator.foreign.*;


public class NativeConfEnergyCalculator implements ConfEnergyCalculator {

	public enum Precision {

		Double(8, MemoryHandles.varHandle(double.class, ByteOrder.nativeOrder())) {

			@Override
			Object fromDouble(double val) {
				return val;
			}

			@Override
			double toDouble(Object val) {
				return (Double)val;
			}
		},

		Single(4, MemoryHandles.varHandle(float.class, ByteOrder.nativeOrder())) {

			@Override
			Object fromDouble(double val) {
				return (float)val;
			}

			@Override
			double toDouble(Object val) {
				return (double)(Float)val;
			}
		};

		public final int bytes;
		public final VarHandle handle;

		Precision(int bytes, VarHandle handle) {
			this.bytes = bytes;
			this.handle = handle;
		}

		abstract Object fromDouble(double val);
		abstract double toDouble(Object val);
	}

	public class Real extends Field {

		private Real() {
			super(precision.bytes);
		}

		public double get() {
			return (double)precision.handle.get(addr);
		}

		public void set(double value) {
			precision.handle.set(addr, precision.fromDouble(value));
		}
	}
	public Real real() {
		return new Real();
	}

	public class RealArray extends Structs.Array {

		public RealArray(long size) {
			super(size, precision.bytes);
		}

		public double get(long i) {
			return (double)precision.handle.get(offset(i));
		}

		public void set(long i, double value) {
			precision.handle.set(offset(i), precision.fromDouble(value));
		}
	}
	public RealArray realarray(long size) {
		return new RealArray(size);
	}


	private static class NativeF64 {

		static {
			Native.register("ConfEcalc_f64");
		}

		// NOTE: there's no way to share these declarations between f64 and f32
		// JNA won't look at methods from base classes

		public static native void print_info();
		public static native void assign(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native double calc(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
	}

	private static class NativeF32 {

		static {
			Native.register("ConfEcalc_f32");
		}

		// NOTE: there's no way to share these declarations between f64 and f32
		// JNA won't look at methods from base classes

		public static native void print_info();
		public static native void assign(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native double calc(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
	}

	private static final long alignment = 8; // bytes

	public final ConfSpace confSpace;
	public final Precision precision;

	private final MemorySegment confSpaceMem;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	static class SConfSpace extends Struct {
		final Uint32 num_pos = uint32();
		final Uint32 max_num_conf_atoms = uint32();
		final Int64 positions_offset = int64();
		final Int64 static_atoms_offset = int64();
	}
	private final SConfSpace confSpaceStruct = new SConfSpace().init(
		"num_pos", "max_num_conf_atoms", "positions_offset", "static_atoms_offset"
	);

	static class SPos extends Struct {
		final Uint32 num_confs = uint32();
		final Uint32 max_num_atoms = uint32();
		final Int64.Array conf_offsets = int64array(Field.UnknownSize);
	}
	private final SPos posStruct = new SPos().init(
		"num_confs", "max_num_atoms", "conf_offsets"
	);

	static class SConf extends Struct {
		final Int64 coords_offset = int64();
	}
	private final SConf confStruct = new SConf().init(
		"coords_offset"
	);

	class SCoords extends Struct {
		final Uint32 num_atoms = uint32();
		final Pad pad = pad(4);
		final RealArray coords = realarray(Field.UnknownSize);
	}
	private final SCoords coordsStruct;

	class SPosInter extends Struct {
		final Uint32 posi1 = uint32();
		final Uint32 posi2 = uint32();
		final Real weight = real();
		final Real offset = real();
	}
	private final SPosInter posInterStruct;

	public NativeConfEnergyCalculator(ConfSpace confSpace, Precision precision) {

		this.confSpace = confSpace;
		this.precision = precision;

		switch (precision) {
			case Single: NativeF32.print_info(); break;
			case Double: NativeF64.print_info(); break;
		}

		// once we know the precision, init the rest of the structs
		coordsStruct = new SCoords().init(
			"num_atoms", "pad", "coords"
		);
		posInterStruct = new SPosInter().init(
			"posi1", "posi2", "weight", "offset"
		);

		Int64.Array posOffsets = int64array(confSpace.positions.length);

		// calculate how much memory we need for the conf space buffer
		// TODO: make this easier to read?
		// TODO: add alignment to bytes(size) for structs with arrays at the end
		// TODO: add alignment to bytes() for fixed arrays
		long bufSize = confSpaceStruct.bytes()
			+ BufWriter.padToAlignment(posOffsets.bytes(), alignment)
			+ Arrays.stream(confSpace.positions)
				.mapToLong(pos ->
					posStruct.bytes(pos.confs.length)
					+ Arrays.stream(pos.confs)
						.mapToLong(conf ->
							confStruct.bytes()
							+ BufWriter.padToAlignment(coordsStruct.bytes(conf.coords.size*3), alignment)
						)
						.sum()
				)
				.sum()
			+ BufWriter.padToAlignment(coordsStruct.bytes(confSpace.staticCoords.size*3), alignment);

		confSpaceMem = MemorySegment.allocateNative(bufSize, alignment);
		BufWriter buf = new BufWriter(confSpaceMem);

		// write the header
		buf.place(confSpaceStruct);
		confSpaceStruct.num_pos.set(confSpace.positions.length);
		confSpaceStruct.max_num_conf_atoms.set(confSpace.maxNumConfAtoms);
		// we'll go back and write the offsets later

		// leave space for the position offsets
		confSpaceStruct.positions_offset.set(buf.pos);
		buf.place(posOffsets);
		buf.skipToAlignment(alignment);

		// write the positions
		for (ConfSpace.Pos pos : confSpace.positions) {
			posOffsets.set(pos.index, buf.pos);
			buf.place(posStruct, pos.confs.length);
			posStruct.num_confs.set(pos.confs.length);
			posStruct.max_num_atoms.set(pos.maxNumAtoms);
			buf.skipToAlignment(alignment);

			// write the confs
			for (ConfSpace.Conf conf : pos.confs) {
				posStruct.conf_offsets.set(conf.index, buf.pos);
				buf.place(confStruct);

				// write the atoms
				confStruct.coords_offset.set(buf.pos);
				buf.place(coordsStruct, conf.coords.size*3);
				coordsStruct.num_atoms.set(conf.coords.size);
				for (int i=0; i<conf.coords.size; i++) {
					coordsStruct.coords.set(i*3,     conf.coords.x(i));
					coordsStruct.coords.set(i*3 + 1, conf.coords.y(i));
					coordsStruct.coords.set(i*3 + 2, conf.coords.z(i));
					// TODO: padding when real_t = float
				}
				buf.skipToAlignment(alignment);
			}
		}

		// write the static atoms
		confSpaceStruct.static_atoms_offset.set(buf.pos);
		buf.place(coordsStruct, confSpace.staticCoords.size*3);
		coordsStruct.num_atoms.set(confSpace.staticCoords.size);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			coordsStruct.coords.set(i*3,     confSpace.staticCoords.x(i));
			coordsStruct.coords.set(i*3 + 1, confSpace.staticCoords.y(i));
			coordsStruct.coords.set(i*3 + 2, confSpace.staticCoords.z(i));
			// TODO: padding when real_t = float
		}
		buf.skipToAlignment(alignment);

		assert(buf.pos == bufSize);

		/* TEMP
		int[] conf = new int[] { 0, 0, 0, 0, 0, 0, 0 };
		var inters = Arrays.asList(
			new PosInter(1, 2, 3.0, 4.0),
			new PosInter(5, 6, 7.0, 8.0)
		);
		try (MemorySegment intersBuf = makeIntersBuf(inters)) {
			NativeF64.calc(
				confSpaceMem.asByteBuffer(),
				conf,
				intersBuf.asByteBuffer(), inters.size()
			);
		}
		*/
	}

	private void assign(int[] conf, ByteBuffer coords) {
		switch (precision) {
			case Single: NativeF32.assign(confSpaceMem.asByteBuffer(), conf, coords); break;
			case Double: NativeF64.assign(confSpaceMem.asByteBuffer(), conf, coords); break;
		}
	}

	public CoordsList assign(int[] conf) {
		try (var coordsMem = MemorySegment.allocateNative(confSpace.maxNumConfAtoms*3*8, alignment)) {

			// run the native code
			assign(conf, coordsMem.asByteBuffer());
			var nums = realarray(confSpace.maxNumConfAtoms);
			nums.setAddress(coordsMem.baseAddress());

			// copy the coords into a CoordsList
			CoordsList coords = new CoordsList(confSpace.maxNumConfAtoms);
			for (int i=0; i<confSpace.maxNumConfAtoms; i++) {
				coords.set(
					i,
					nums.get(i*3),
					nums.get(i*3 + 1),
					nums.get(i*3 + 2)
				);
			}
			return coords;
		}
	}

	private MemorySegment makeIntersBuf(List<PosInter> inters) {
		MemorySegment mem = MemorySegment.allocateNative(posInterStruct.bytes()*inters.size(), alignment);
		BufWriter buf = new BufWriter(mem);
		for (var inter : inters) {
			buf.place(posInterStruct);
			posInterStruct.posi1.set(inter.posi1);
			posInterStruct.posi2.set(inter.posi2);
			posInterStruct.weight.set(inter.weight);
			posInterStruct.offset.set(inter.offset);
		}
		return mem;
	}

	@Override
	public void close() {
		confSpaceMem.close();
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {
		// TODO
		throw new Error();
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {
		// TODO
		throw new Error();
	}
}
