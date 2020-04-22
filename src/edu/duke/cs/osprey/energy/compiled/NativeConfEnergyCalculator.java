package edu.duke.cs.osprey.energy.compiled;

import com.sun.jna.*;

import static edu.duke.cs.osprey.gpu.Structs.*;

import java.lang.invoke.VarHandle;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.List;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.gpu.BufWriter;
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


	private static class NativeLib {

		static {
			Native.register("ConfEcalc");
		}

		public static native int version_major();
		public static native int version_minor();
		public static native void assign_f32(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native void assign_f64(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		//public static native double calc(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
	}

	private static final long alignment = 8; // bytes

	public final ConfSpace confSpace;
	public final Precision precision;

	private final MemorySegment confSpaceMem;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	static class SConfSpace extends Struct {
		final Int32 num_pos = int32();
		final Int32 max_num_conf_atoms = int32();
		final Int64 positions_offset = int64();
		final Int64 static_atoms_offset = int64();
	}
	private final SConfSpace confSpaceStruct = new SConfSpace().init(
		24, "num_pos", "max_num_conf_atoms", "positions_offset", "static_atoms_offset"
	);

	static class SPos extends Struct {
		final Int32 num_confs = int32();
		final Int32 max_num_atoms = int32();
	}
	private final SPos posStruct = new SPos().init(
		8, "num_confs", "max_num_atoms"
	);

	static class SConf extends Struct {
		final Int64 atoms_offset = int64();
	}
	private final SConf confStruct = new SConf().init(
		8, "atoms_offset"
	);

	static class SAtoms extends Struct {
		final Int32 num_atoms = int32();
		final Pad pad = pad(4);
		final Int64 coords = int64();
	}
	private final SAtoms atomsStruct = new SAtoms().init(
		16, "num_atoms", "pad", "coords"
	);

	class SReal3 extends Struct {
		final Real x = real();
		final Real y = real();
		final Real z = real();
		final Pad pad = pad(switch (precision) {
			case Single -> 4;
			case Double -> 0;
		});
	}
	private final SReal3 real3Struct;

	class SPosInter extends Struct {
		final Int32 posi1 = int32();
		final Int32 posi2 = int32();
		final Real weight = real();
		final Real offset = real();
	}
	private final SPosInter posInterStruct;

	public NativeConfEnergyCalculator(ConfSpace confSpace, Precision precision) {

		this.confSpace = confSpace;
		this.precision = precision;

		// once we know the precision, init the rest of the structs
		real3Struct = new SReal3().init(
			switch (precision) {
				case Single -> 16;
				case Double -> 24;
			}, "x", "y", "z", "pad"
		);
		posInterStruct = new SPosInter().init(
			switch (precision) {
				case Single -> 16;
				case Double -> 24;
			}, "posi1", "posi2", "weight", "offset"
		);

		Int64.Array posOffsets = int64array();
		Int64.Array confOffsets = int64array();

		// calculate how much memory we need for the conf space buffer
		long bufSize = confSpaceStruct.bytes()
			+ posOffsets.bytes(confSpace.positions.length)
			+ sum(confSpace.positions, pos ->
				posStruct.bytes()
				+ confOffsets.bytes(pos.confs.length)
				+ sum(pos.confs, conf ->
					confStruct.bytes()
					+ atomsStruct.bytes()
					+ real3Struct.bytes()*conf.coords.size
				)
			)
			+ atomsStruct.bytes()
			+ real3Struct.bytes()*confSpace.staticCoords.size;

		// allocate the buffer for the conf space
		confSpaceMem = MemorySegment.allocateNative(bufSize, alignment);
		BufWriter buf = new BufWriter(confSpaceMem);

		// write the header
		buf.place(confSpaceStruct);
		confSpaceStruct.num_pos.set(confSpace.positions.length);
		confSpaceStruct.max_num_conf_atoms.set(confSpace.maxNumConfAtoms);
		// we'll go back and write the offsets later

		// leave space for the position offsets
		confSpaceStruct.positions_offset.set(buf.pos);
		buf.place(posOffsets, confSpace.positions.length);

		// write the positions
		for (ConfSpace.Pos pos : confSpace.positions) {
			posOffsets.set(pos.index, buf.pos);
			buf.place(posStruct);
			posStruct.num_confs.set(pos.confs.length);
			posStruct.max_num_atoms.set(pos.maxNumAtoms);

			// put the conf offsets
			buf.place(confOffsets, pos.confs.length);
			// we'll go back and write them later though

			// write the confs
			for (ConfSpace.Conf conf : pos.confs) {
				confOffsets.set(conf.index, buf.pos);
				buf.place(confStruct);

				// write the atoms
				confStruct.atoms_offset.set(buf.pos);
				buf.place(atomsStruct);
				atomsStruct.num_atoms.set(conf.coords.size);
				atomsStruct.coords.set(0);
				for (int i=0; i<conf.coords.size; i++) {
					buf.place(real3Struct);
					real3Struct.x.set(conf.coords.x(i));
					real3Struct.y.set(conf.coords.y(i));
					real3Struct.z.set(conf.coords.z(i));
				}
			}
		}

		// write the static atoms
		confSpaceStruct.static_atoms_offset.set(buf.pos);
		buf.place(atomsStruct);
		atomsStruct.num_atoms.set(confSpace.staticCoords.size);
		atomsStruct.coords.set(0);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			buf.place(real3Struct);
			real3Struct.x.set(confSpace.staticCoords.x(i));
			real3Struct.y.set(confSpace.staticCoords.y(i));
			real3Struct.z.set(confSpace.staticCoords.z(i));
		}

		assert(buf.pos == bufSize);

		/* TODO: position interactions
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

	public String version() {
		return String.format("%d.%d", NativeLib.version_major(), NativeLib.version_minor());
	}

	private void assign(int[] conf, ByteBuffer coords) {
		switch (precision) {
			case Single: NativeLib.assign_f32(confSpaceMem.asByteBuffer(), conf, coords); break;
			case Double: NativeLib.assign_f64(confSpaceMem.asByteBuffer(), conf, coords); break;
		}
	}

	public CoordsList assign(int[] conf) {
		try (var coordsMem = MemorySegment.allocateNative(real3Struct.bytes()*confSpace.maxNumConfAtoms, alignment)) {

			// run the native code
			assign(conf, coordsMem.asByteBuffer());

			// copy the coords into a CoordsList
			CoordsList coords = new CoordsList(confSpace.maxNumConfAtoms);
			var addr = coordsMem.baseAddress();
			for (int i=0; i<confSpace.maxNumConfAtoms; i++) {
				real3Struct.setAddress(addr.addOffset(i*real3Struct.bytes()));
				coords.set(
					i,
					real3Struct.x.get(),
					real3Struct.y.get(),
					real3Struct.z.get()
				);
			}
			return coords;
		}
	}

	// TODO: update me!
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
