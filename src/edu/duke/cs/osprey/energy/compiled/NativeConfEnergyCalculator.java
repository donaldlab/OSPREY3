package edu.duke.cs.osprey.energy.compiled;

import com.sun.jna.*;

import static edu.duke.cs.osprey.gpu.Structs.*;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.CoordsList;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.gpu.BufWriter;
import jdk.incubator.foreign.*;


public class NativeConfEnergyCalculator implements ConfEnergyCalculator {

	public enum Precision {

		Double(8),
		Float(4);

		public final int bytes;

		Precision(int bytes) {
			this.bytes = bytes;
		}
	}

	private static class NativeF64 {

		static {
			Native.register("ConfEcalc_f64");
		}

		public static native int version_major();
		public static native int version_minor();
		public static native void assign(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native double calc(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
	}

	private static class NativeF32 {

		static {
			Native.register("ConfEcalc_f32");
		}

		// TODO
	}

	private static final long alignment = 8; // bytes

	public final ConfSpace confSpace;
	public final Precision precision;

	private final MemorySegment confSpaceMem;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	public static class SConfSpace extends Struct {
		public final Uint32 num_pos = uint32();
		public final Uint32 max_num_conf_atoms = uint32();
		public final Int64 positions_offset = int64();
		public final Int64 static_atoms_offset = int64();
	}
	private final SConfSpace confSpaceStruct = new SConfSpace().init(
		"num_pos", "max_num_conf_atoms", "positions_offset", "static_atoms_offset"
	);

	public static class SPos extends Struct {
		public final Uint32 num_confs = uint32();
		public final Uint32 max_num_atoms = uint32();
		public final Int64.Array conf_offsets = int64array(Field.UnknownSize);
	}
	private final SPos posStruct = new SPos().init(
		"num_confs", "max_num_atoms", "conf_offsets"
	);

	public static class SConf extends Struct {
		public final Int64 coords_offset = int64();
	}
	private final SConf confStruct = new SConf().init(
		"coords_offset"
	);

	public static class SCoords extends Struct {
		public final Uint32 num_atoms = uint32();
		public final Pad pad = pad(4);
		public final Float64.Array coords = float64array(Field.UnknownSize);
	}
	private final SCoords coordsStruct = new SCoords().init(
		"num_atoms", "pad", "coords"
	);

	public static class SPosInter extends Struct {
		public final Uint32 posi1 = uint32();
		public final Uint32 posi2 = uint32();
		public final Float64 weight = float64();
		public final Float64 offset = float64();
	}
	private final SPosInter posInterStruct = new SPosInter().init(
		"posi1", "posi2", "weight", "offset"
	);

	public NativeConfEnergyCalculator(ConfSpace confSpace, Precision precision) {

		this.confSpace = confSpace;
		this.precision = precision;

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

		// TEMP
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
	}

	public CoordsList assign(int[] conf) {
		try (var coordsMem = MemorySegment.allocateNative(confSpace.maxNumConfAtoms*3*8, alignment)) {

			// run the native code
			NativeF64.assign(confSpaceMem.asByteBuffer(), conf, coordsMem.asByteBuffer());
			Float64.Array nums = float64array(confSpace.maxNumConfAtoms);
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
