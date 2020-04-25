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
		public static native float calc_energy_amber_eef1_f32(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
		public static native double calc_energy_amber_eef1_f64(ByteBuffer conf_space, int[] conf, ByteBuffer inters, int inters_size);
	}

	private interface ForcefieldsImpl {
		double calcEnergy(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, int intersSize);
		long paramsBytes();
		void writeParams(ConfSpace confSpace, BufWriter buf);
		long staticStaticBytes(ConfSpace confSpace);
		long staticPosBytes(ConfSpace confSpace, int posi1, int fragi1);
		long posBytes(ConfSpace confSpace, int posi1, int fragi1);
		long posPosBytes(ConfSpace confSpace, int posi1, int fragi1, int posi2, int fragi2);
		void writeStaticStatic(ConfSpace confSpace, BufWriter buf);
		void writeStaticPos(ConfSpace confSpace, int posi1, int fragi1, BufWriter buf);
		void writePos(ConfSpace confSpace, int posi1, int fragi1, BufWriter buf);
		void writePosPos(ConfSpace confSpace, int posi1, int fragi1, int posi2, int fragi2, BufWriter buf);
	}

	private class AmberEef1 implements ForcefieldsImpl {

		class SParamsAmberEef1 extends Struct {
			Bool distance_dependent_dielectric = bool();
			Pad pad = pad(7);
		}
		SParamsAmberEef1 paramsStruct = new SParamsAmberEef1().init(
			8, "distance_dependent_dielectric", "pad"
		);

		class SAtomPairs extends Struct {
			final Int32 num_amber = int32();
			final Int32 num_eef1 = int32();
		}
		SAtomPairs atomPairsStruct = new SAtomPairs().init(
			8, "num_amber", "num_eef1"
		);

		class SAtomPairAmber extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			final Real esQ = real();
			final Real vdwA = real();
			final Real vdwB = real();
			final Pad pad = pad(switch (precision) {
				case Single -> 4;
				case Double -> 0;
			});

			void setParams(double[] params) {
				esQ.set(params[0]);
				vdwA.set(params[1]);
				vdwB.set(params[2]);
			}
		}
		SAtomPairAmber amberStruct;

		class SAtomPairEef1 extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			final Real vdwRadius1 = real();
			final Real lambda1 = real();
			final Real vdwRadius2 = real();
			final Real lambda2 = real();
			final Real alpha1 = real();
			final Real alpha2 = real();

			void setParams(double[] params) {
				vdwRadius1.set(params[0]);
				lambda1.set(params[1]);
				vdwRadius2.set(params[2]);
				lambda2.set(params[3]);
				alpha1.set(params[4]);
				alpha2.set(params[5]);
			}
		}
		SAtomPairEef1 eef1Struct;

		AmberEef1(Precision precision) {

			// once we know the precision, init the rest of the structs
			amberStruct = new SAtomPairAmber().init(
				switch (precision) {
					case Single -> 24;
					case Double -> 32;
				}, "atomi1", "atomi2", "esQ", "vdwA", "vdwB", "pad"
			);

			eef1Struct = new SAtomPairEef1().init(
				switch (precision) {
					case Single -> 32;
					case Double -> 56;
				}, "atomi1", "atomi2", "vdwRadius1", "lambda1", "vdwRadius2", "lambda2", "alpha1", "alpha2"
			);
		}

		@Override
		public double calcEnergy(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, int intersSize) {
			return switch (precision) {
				case Single -> NativeLib.calc_energy_amber_eef1_f32(confSpaceBuf, conf, intersBuf, intersSize);
				case Double -> NativeLib.calc_energy_amber_eef1_f64(confSpaceBuf, conf, intersBuf, intersSize);
			};
		}

		@Override
		public long paramsBytes() {
			return paramsStruct.bytes();
		}

		@Override
		public void writeParams(ConfSpace confSpace, BufWriter buf) {

			buf.place(paramsStruct);

			// write the amber params
			AmberEnergyCalculator.Settings amberSettings = ((AmberEnergyCalculator)confSpace.ecalcs[0]).settings;
			paramsStruct.distance_dependent_dielectric.set(amberSettings.distanceDependentDielectric);

			// no EEF1 params to write
		}

		@Override
		public long staticStaticBytes(ConfSpace confSpace) {
			return atomPairsStruct.bytes()
				+ confSpace.indicesStatic(0).size()*amberStruct.bytes()
				+ confSpace.indicesStatic(1).size()*eef1Struct.bytes();
		}

		@Override
		public long staticPosBytes(ConfSpace confSpace, int posi1, int fragi1) {
			return atomPairsStruct.bytes()
				+ confSpace.indicesSinglesByFrag(0, posi1, fragi1).sizeStatics()*amberStruct.bytes()
				+ confSpace.indicesSinglesByFrag(1, posi1, fragi1).sizeStatics()*eef1Struct.bytes();
		}

		@Override
		public long posBytes(ConfSpace confSpace, int posi1, int fragi1) {
			return atomPairsStruct.bytes()
				+ confSpace.indicesSinglesByFrag(0, posi1, fragi1).sizeInternals()*amberStruct.bytes()
				+ confSpace.indicesSinglesByFrag(1, posi1, fragi1).sizeInternals()*eef1Struct.bytes();
		}

		@Override
		public long posPosBytes(ConfSpace confSpace, int posi1, int fragi1, int posi2, int fragi2) {
			return atomPairsStruct.bytes()
				+ confSpace.indicesPairsByFrags(0, posi1, fragi1, posi2, fragi2).size()*amberStruct.bytes()
				+ confSpace.indicesPairsByFrags(1, posi1, fragi1, posi2, fragi2).size()*eef1Struct.bytes();
		}

		@Override
		public void writeStaticStatic(ConfSpace confSpace, BufWriter buf) {

			ConfSpace.IndicesStatic amberIndices = confSpace.indicesStatic(0);
			ConfSpace.IndicesStatic eef1Indices = confSpace.indicesStatic(1);

			buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(amberIndices.size());
			atomPairsStruct.num_eef1.set(eef1Indices.size());

			for (int i=0; i<amberIndices.size(); i++) {
				buf.place(amberStruct);
				amberStruct.atomi1.set(confSpace.getStaticAtomIndex(amberIndices.getStaticAtom1Index(i)));
				amberStruct.atomi2.set(confSpace.getStaticAtomIndex(amberIndices.getStaticAtom2Index(i)));
				amberStruct.setParams(confSpace.ffparams(0, amberIndices.getParamsIndex(i)));
			}

			for (int i=0; i<eef1Indices.size(); i++) {
				buf.place(eef1Struct);
				eef1Struct.atomi1.set(confSpace.getStaticAtomIndex(eef1Indices.getStaticAtom1Index(i)));
				eef1Struct.atomi2.set(confSpace.getStaticAtomIndex(eef1Indices.getStaticAtom2Index(i)));
				eef1Struct.setParams(confSpace.ffparams(1, eef1Indices.getParamsIndex(i)));
			}
		}

		@Override
		public void writeStaticPos(ConfSpace confSpace, int posi1, int fragi1, BufWriter buf) {

			ConfSpace.IndicesSingle amberIndices = confSpace.indicesSinglesByFrag(0, posi1, fragi1);
			ConfSpace.IndicesSingle eef1Indices = confSpace.indicesSinglesByFrag(1, posi1, fragi1);

			buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(amberIndices.sizeStatics());
			atomPairsStruct.num_eef1.set(eef1Indices.sizeStatics());

			for (int i=0; i<amberIndices.sizeStatics(); i++) {
				buf.place(amberStruct);
				amberStruct.atomi1.set(confSpace.getStaticAtomIndex(amberIndices.getStaticStaticAtomIndex(i)));
				amberStruct.atomi2.set(confSpace.getConfAtomIndex(posi1, amberIndices.getStaticConfAtomIndex(i)));
				amberStruct.setParams(confSpace.ffparams(0, amberIndices.getStaticParamsIndex(i)));
			}

			for (int i=0; i<eef1Indices.sizeStatics(); i++) {
				buf.place(eef1Struct);
				eef1Struct.atomi1.set(confSpace.getStaticAtomIndex(eef1Indices.getStaticStaticAtomIndex(i)));
				eef1Struct.atomi2.set(confSpace.getConfAtomIndex(posi1, eef1Indices.getStaticConfAtomIndex(i)));
				eef1Struct.setParams(confSpace.ffparams(1, eef1Indices.getStaticParamsIndex(i)));
			}
		}

		@Override
		public void writePos(ConfSpace confSpace, int posi1, int fragi1, BufWriter buf) {

			ConfSpace.IndicesSingle amberIndices = confSpace.indicesSinglesByFrag(0, posi1, fragi1);
			ConfSpace.IndicesSingle eef1Indices = confSpace.indicesSinglesByFrag(1, posi1, fragi1);

			buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(amberIndices.sizeInternals());
			atomPairsStruct.num_eef1.set(eef1Indices.sizeInternals());

			for (int i=0; i<amberIndices.sizeInternals(); i++) {
				buf.place(amberStruct);
				amberStruct.atomi1.set(confSpace.getConfAtomIndex(posi1, amberIndices.getInternalConfAtom1Index(i)));
				amberStruct.atomi2.set(confSpace.getConfAtomIndex(posi1, amberIndices.getInternalConfAtom2Index(i)));
				amberStruct.setParams(confSpace.ffparams(0, amberIndices.getInternalParamsIndex(i)));
			}

			for (int i=0; i<eef1Indices.sizeInternals(); i++) {
				buf.place(eef1Struct);
				eef1Struct.atomi1.set(confSpace.getConfAtomIndex(posi1, eef1Indices.getInternalConfAtom1Index(i)));
				eef1Struct.atomi2.set(confSpace.getConfAtomIndex(posi1, eef1Indices.getInternalConfAtom2Index(i)));
				eef1Struct.setParams(confSpace.ffparams(1, eef1Indices.getInternalParamsIndex(i)));
			}
		}

		@Override
		public void writePosPos(ConfSpace confSpace, int posi1, int fragi1, int posi2, int fragi2, BufWriter buf) {

			ConfSpace.IndicesPair amberIndices = confSpace.indicesPairsByFrags(0, posi1, fragi1, posi2, fragi2);
			ConfSpace.IndicesPair eef1Indices = confSpace.indicesPairsByFrags(1, posi1, fragi1, posi2, fragi2);

			buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(amberIndices.size());
			atomPairsStruct.num_eef1.set(eef1Indices.size());

			for (int i=0; i<amberIndices.size(); i++) {
				buf.place(amberStruct);
				amberStruct.atomi1.set(confSpace.getConfAtomIndex(posi1, amberIndices.getConfAtom1Index(i)));
				amberStruct.atomi2.set(confSpace.getConfAtomIndex(posi2, amberIndices.getConfAtom2Index(i)));
				amberStruct.setParams(confSpace.ffparams(0, amberIndices.getParamsIndex(i)));
			}

			for (int i=0; i<eef1Indices.size(); i++) {
				buf.place(eef1Struct);
				eef1Struct.atomi1.set(confSpace.getConfAtomIndex(posi1, eef1Indices.getConfAtom1Index(i)));
				eef1Struct.atomi2.set(confSpace.getConfAtomIndex(posi2, eef1Indices.getConfAtom2Index(i)));
				eef1Struct.setParams(confSpace.ffparams(1, eef1Indices.getParamsIndex(i)));
			}
		}
	}

	public final ConfSpace confSpace;
	public final Precision precision;
	public final ForcefieldsImpl forcefieldsImpl;

	private final MemorySegment confSpaceMem;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	class SConfSpace extends Struct {
		final Int32 num_pos = int32();
		final Int32 max_num_conf_atoms = int32();
		final Int64 positions_offset = int64();
		final Int64 static_atoms_offset = int64();
		final Int64 params_offset = int64();
		final Int64 pos_pairs_offset = int64();
		final Real static_energy = real();
		final Pad pad = pad(switch (precision) {
			case Single -> 4;
			case Double -> 0;
		});
	}
	private final SConfSpace confSpaceStruct;

	static class SPos extends Struct {
		final Int32 num_confs = int32();
		final Int32 max_num_atoms = int32();
		final Int32 num_frags = int32();
		final Pad pad = pad(4);
	}
	private final SPos posStruct = new SPos().init(
		16, "num_confs", "max_num_atoms", "num_frags", "pad"
	);

	class SConf extends Struct {
		final Int64 atoms_offset = int64();
		final Int32 frag_index = int32();
		final Pad pad = pad(switch (precision) {
			case Single -> 0;
			case Double -> 4;
		});
		final Real internal_energy = real();
	}
	private final SConf confStruct;

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

		// find the forcefield implementation, or die trying
		EnergyCalculator.Type[] ecalcTypes = Arrays.stream(confSpace.ecalcs)
			.map(EnergyCalculator::type)
			.toArray(EnergyCalculator.Type[]::new);
		if (ecalcTypes.length == 2 && ecalcTypes[0] == AmberEnergyCalculator.type && ecalcTypes[1] == EEF1EnergyCalculator.type) {
			forcefieldsImpl = new AmberEef1(precision);
		} else {
			throw new IllegalArgumentException("No native implementation for forcefields: " + Arrays.toString(ecalcTypes));
		}

		// once we know the precision, init the rest of the structs
		confSpaceStruct = new SConfSpace().init(
			switch (precision) {
				case Single -> 48;
				case Double -> 48;
			},
			"num_pos", "max_num_conf_atoms",
			"positions_offset", "static_atoms_offset", "params_offset", "pos_pairs_offset",
			"static_energy", "pad"
		);
		confStruct = new SConf().init(
			switch (precision) {
				case Single -> 16;
				case Double -> 24;
			}, "atoms_offset", "frag_index", "pad", "internal_energy"
		);
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
		Int64.Array posPairOffsets = int64array();

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
			+ real3Struct.bytes()*confSpace.staticCoords.size
			+ forcefieldsImpl.paramsBytes();

		// add space for the pos pair offsets
		int numPosPairs =
			1 // static-static
			+ confSpace.numPos() // static-pos
			+ confSpace.numPos() // pos
			+ confSpace.numPos()*(confSpace.numPos() - 1)/2; // pos-pos
		bufSize += posPairOffsets.bytes(numPosPairs);

		// add space for the static-static pairs
		bufSize += forcefieldsImpl.staticStaticBytes(confSpace);

		// add space for the static-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			bufSize += Int64.bytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.staticPosBytes(confSpace, posi1, fragi1);
			}
		}

		// add space for the pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			bufSize += Int64.bytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.posBytes(confSpace, posi1, fragi1);

			}
		}

		// add space for the pos-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				bufSize += Int64.bytes*confSpace.numFrag(posi1)*confSpace.numFrag(posi2);
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						bufSize += forcefieldsImpl.posPosBytes(confSpace, posi1, fragi1, posi2, fragi2);
					}
				}
			}
		}

		// allocate the buffer for the conf space
		confSpaceMem = MemorySegment.allocateNative(bufSize);
		BufWriter buf = new BufWriter(confSpaceMem);

		// write the header
		buf.place(confSpaceStruct);
		confSpaceStruct.num_pos.set(confSpace.positions.length);
		confSpaceStruct.max_num_conf_atoms.set(confSpace.maxNumConfAtoms);
		// we'll go back and write the offsets later
		confSpaceStruct.static_energy.set(Arrays.stream(confSpace.staticEnergies).sum());

		// leave space for the position offsets
		confSpaceStruct.positions_offset.set(buf.pos);
		buf.place(posOffsets, confSpace.positions.length);

		// write the positions
		for (ConfSpace.Pos pos : confSpace.positions) {
			posOffsets.set(pos.index, buf.pos);
			buf.place(posStruct);
			posStruct.num_confs.set(pos.confs.length);
			posStruct.max_num_atoms.set(pos.maxNumAtoms);
			posStruct.num_frags.set(pos.numFrags);

			// put the conf offsets
			buf.place(confOffsets, pos.confs.length);
			// we'll go back and write them later though

			// write the confs
			for (ConfSpace.Conf conf : pos.confs) {
				confOffsets.set(conf.index, buf.pos);
				buf.place(confStruct);
				confStruct.frag_index.set(conf.fragIndex);
				confStruct.internal_energy.set(Arrays.stream(conf.energies).sum());

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

		// write the forcefield params
		confSpaceStruct.params_offset.set(buf.pos);
		forcefieldsImpl.writeParams(confSpace, buf);

		// write the pos pairs
		confSpaceStruct.pos_pairs_offset.set(buf.pos);
		buf.place(posPairOffsets, numPosPairs);

		// write the static-static pair
		posPairOffsets.set(0, buf.pos);
		forcefieldsImpl.writeStaticStatic(confSpace, buf);

		// write the static-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(1 + posi1, buf.pos);
			Int64.Array fragOffsets = new Int64.Array();
			buf.place(fragOffsets, confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragi1, buf.pos);
				forcefieldsImpl.writeStaticPos(confSpace, posi1, fragi1, buf);
			}
		}

		// write the pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(1 + confSpace.positions.length + posi1, buf.pos);
			Int64.Array fragOffsets = new Int64.Array();
			buf.place(fragOffsets, confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragi1, buf.pos);
				forcefieldsImpl.writePos(confSpace, posi1, fragi1, buf);
			}
		}

		// write the pos-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				posPairOffsets.set(1 + 2*confSpace.positions.length + posi1*(posi1 - 1)/2 + posi2, buf.pos);
				Int64.Array fragOffsets = new Int64.Array();
				buf.place(fragOffsets, confSpace.numFrag(posi1)*confSpace.numFrag(posi2));
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						fragOffsets.set(fragi1*confSpace.numFrag(posi2) + fragi2, buf.pos);
						forcefieldsImpl.writePosPos(confSpace, posi1, fragi1, posi2, fragi2, buf);
					}
				}
			}
		}

		assert(buf.pos == bufSize);
	}

	public String version() {
		return String.format("%d.%d", NativeLib.version_major(), NativeLib.version_minor());
	}

	private void assign(int[] conf, ByteBuffer coords) {
		switch (precision) {
			case Single -> NativeLib.assign_f32(confSpaceMem.asByteBuffer(), conf, coords);
			case Double -> NativeLib.assign_f64(confSpaceMem.asByteBuffer(), conf, coords);
		}
	}

	public CoordsList assign(int[] conf) {
		try (var coordsMem = MemorySegment.allocateNative(real3Struct.bytes()*confSpace.maxNumConfAtoms)) {

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
	public double calcEnergy(int[] conf, List<PosInter> inters) {
		try (MemorySegment intersMem = makeIntersMem(inters)) {
			return forcefieldsImpl.calcEnergy(confSpaceMem.asByteBuffer(), conf, intersMem.asByteBuffer(), inters.size());
		}
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {
		// TODO
		throw new Error();
	}

	private MemorySegment makeIntersMem(List<PosInter> inters) {
		MemorySegment mem = MemorySegment.allocateNative(posInterStruct.bytes()*inters.size());
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
}
