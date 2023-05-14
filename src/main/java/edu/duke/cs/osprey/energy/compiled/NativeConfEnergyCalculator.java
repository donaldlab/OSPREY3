package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.sun.jna.*;

import static edu.duke.cs.osprey.gpu.Structs.*;

import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.confspace.compiled.motions.TranslationRotation;
import edu.duke.cs.osprey.gpu.MemoryBuffer;
import edu.duke.cs.osprey.gpu.Precision;


/**
 * A reference implementation of a conformation minimizer and energy calculator in C++, rather than Java.
 *
 * It's not terribly optimized, but it's already faster than the original implementation in Java.
 */
/*
 * This is the Java side, which is mostly responsible for converting the Java objects into
 * something C++ can understand (and dealing with memory layouts, etc), calling the C++ library
 * over Java's FFI, and then convertering the C++ results back into something Java can understand.
 *
 * The FFI in this case is based on JNA. As java.lang.foreign is (currently) in preview in JDK 19, a MemoryBuffer
 * is used to construct the buffer passed to the native code.
 *
 * This implementation is correct (as best as I can tell), but it's very naive and totally unoptimized.
 */
public class NativeConfEnergyCalculator implements ConfEnergyCalculator {

	private static class NativeLib {

		static {
			Native.register("ConfEcalc");
		}

		public static native int version_major();
		public static native int version_minor();
		public static native void assign_f32(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native void assign_f64(ByteBuffer conf_space, int[] conf, ByteBuffer out);
		public static native float calc_amber_eef1_f32(ByteBuffer conf_space, int[] conf, ByteBuffer inters, ByteBuffer out_coords);
		public static native double calc_amber_eef1_f64(ByteBuffer conf_space, int[] conf, ByteBuffer inters, ByteBuffer out_coords);
		public static native float minimize_amber_eef1_f32(ByteBuffer conf_space, int[] conf, ByteBuffer inters, ByteBuffer out_coords, ByteBuffer out_dofs);
		public static native double minimize_amber_eef1_f64(ByteBuffer conf_space, int[] conf, ByteBuffer inters, ByteBuffer out_coords, ByteBuffer out_dofs);
	}

	private interface ForcefieldsImpl {
		double calc(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords);
		double minimize(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords, ByteBuffer dofs);
		long paramsBytes();
		void writeParams(MemoryBuffer buf);
		long staticStaticBytes();
		long staticPosBytes(int posi1, int fragi1);
		long posBytes(int posi1, int fragi1);
		long posPosBytes(int posi1, int fragi1, int posi2, int fragi2);
		void writeStaticStatic(MemoryBuffer buf);
		void writeStaticPos(int posi1, int fragi1, MemoryBuffer buf);
		void writePos(int posi1, int fragi1, MemoryBuffer buf);
		void writePosPos(int posi1, int fragi1, int posi2, int fragi2, MemoryBuffer buf);
	}

	private interface AtomPairWriter {
		int size();
		int atomi1(int i);
		int atomi2(int i);
		double[] params(int i);
	}

	private class AmberEef1 implements ForcefieldsImpl {

		class SParamsAmberEef1 extends Struct {
			Bool distance_dependent_dielectric = bool();
			Pad pad = pad(7);
			void init() {
				init(8, "distance_dependent_dielectric", "pad");
			}
		}
		final SParamsAmberEef1 paramsStruct = new SParamsAmberEef1();

		class SAtomPairs extends Struct {
			final Int32 num_amber = int32();
			final Int32 num_eef1 = int32();
			void init() {
				init(8, "num_amber", "num_eef1");
			}
		}
		final SAtomPairs atomPairsStruct = new SAtomPairs();

		class SAtomPairAmber extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			final Real esQ = real(precision);
			final Real vdwA = real(precision);
			final Real vdwB = real(precision);
			final Pad pad = pad(precision.map(4, 0));

			void init() {
				init(
					precision.map(24, 32),
				"atomi1", "atomi2", "esQ", "vdwA", "vdwB", "pad"
				);
			}

			void setParams(MemoryBuffer buf, double[] params) {
				esQ.set(buf, params[0]);
				vdwA.set(buf, params[1]);
				vdwB.set(buf, params[2]);
			}
		}
		final SAtomPairAmber amberStruct = new SAtomPairAmber();

		class SAtomPairEef1 extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			final Real vdwRadius1 = real(precision);
			final Real lambda1 = real(precision);
			final Real vdwRadius2 = real(precision);
			final Real lambda2 = real(precision);
			final Real alpha1 = real(precision);
			final Real alpha2 = real(precision);

			void init() {
				init(
					precision.map(32, 56),
					"atomi1", "atomi2", "vdwRadius1", "lambda1", "vdwRadius2", "lambda2", "alpha1", "alpha2"
				);
			}

			void setParams(MemoryBuffer buf, double[] params) {
				vdwRadius1.set(buf, params[0]);
				lambda1.set(buf, params[1]);
				vdwRadius2.set(buf, params[2]);
				lambda2.set(buf, params[3]);
				alpha1.set(buf, params[4]);
				alpha2.set(buf, params[5]);
			}
		}
		final SAtomPairEef1 eef1Struct = new SAtomPairEef1();

		AmberEef1() {

			// once we know the precision, init the structs
			paramsStruct.init();
			atomPairsStruct.init();
			amberStruct.init();
			eef1Struct.init();
		}

		@Override
		public double calc(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords) {
			return switch (precision) {
				case Float32 -> NativeLib.calc_amber_eef1_f32(confSpaceBuf, conf, intersBuf, coords);
				case Float64 -> NativeLib.calc_amber_eef1_f64(confSpaceBuf, conf, intersBuf, coords);
			};
		}

		@Override
		public double minimize(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords, ByteBuffer dofs) {
			return switch (precision) {
				case Float32 -> NativeLib.minimize_amber_eef1_f32(confSpaceBuf, conf, intersBuf, coords, dofs);
				case Float64 -> NativeLib.minimize_amber_eef1_f64(confSpaceBuf, conf, intersBuf, coords, dofs);
			};
		}

		@Override
		public long paramsBytes() {
			return paramsStruct.bytes();
		}

		@Override
		public void writeParams(MemoryBuffer buf) {

			var addr = buf.place(paramsStruct);

			// write the amber params
			AmberEnergyCalculator.Settings amberSettings = ((AmberEnergyCalculator)confSpace.ecalcs[0]).settings;
			paramsStruct.distance_dependent_dielectric.set(addr, amberSettings.distanceDependentDielectric);

			// no EEF1 params to write
		}

		private long atomPairsBytes(int numAmber, int numEef1) {
			return atomPairsStruct.bytes()
				+ numAmber*amberStruct.bytes()
				+ numEef1*eef1Struct.bytes();
		}

		@Override
		public long staticStaticBytes() {
			return atomPairsBytes(
				confSpace.indicesStatic(0).size(),
				confSpace.indicesStatic(1).size()
			);
		}

		@Override
		public long staticPosBytes(int posi1, int fragi1) {
			return atomPairsBytes(
				confSpace.indicesSinglesByFrag(0, posi1, fragi1).sizeStatics(),
				confSpace.indicesSinglesByFrag(1, posi1, fragi1).sizeStatics()
			);
		}

		@Override
		public long posBytes(int posi1, int fragi1) {
			return atomPairsBytes(
				confSpace.indicesSinglesByFrag(0, posi1, fragi1).sizeInternals(),
				confSpace.indicesSinglesByFrag(1, posi1, fragi1).sizeInternals()
			);
		}

		@Override
		public long posPosBytes(int posi1, int fragi1, int posi2, int fragi2) {
			return atomPairsBytes(
				confSpace.indicesPairsByFrags(0, posi1, fragi1, posi2, fragi2).size(),
				confSpace.indicesPairsByFrags(1, posi1, fragi1, posi2, fragi2).size()
			);
		}

		private void writeAtomPairs(AtomPairWriter amber, AtomPairWriter eef1, MemoryBuffer buf) {

			long firstPos = buf.getPos();

			var atomPairsAddr = buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(atomPairsAddr, amber.size());
			atomPairsStruct.num_eef1.set(atomPairsAddr, eef1.size());

			for (int i=0; i<amber.size(); i++) {
				var addr = buf.place(amberStruct);
				int atomi1 = amber.atomi1(i);
				int atomi2 = amber.atomi2(i);
				assert (atomi1 != atomi2);
				amberStruct.atomi1.set(addr, atomi1);
				amberStruct.atomi2.set(addr, atomi2);
				amberStruct.setParams(addr, amber.params(i));
			}

			for (int i=0; i<eef1.size(); i++) {
				var addr = buf.place(eef1Struct);
				int atomi1 = eef1.atomi1(i);
				int atomi2 = eef1.atomi2(i);
				assert (atomi1 != atomi2);
				eef1Struct.atomi1.set(addr, atomi1);
				eef1Struct.atomi2.set(addr, atomi2);
				eef1Struct.setParams(addr, eef1.params(i));
			}

			assert (buf.getPos() - firstPos == atomPairsBytes(amber.size(), eef1.size()))
				: String.format("overshot by %d bytes", buf.getPos() - firstPos - atomPairsBytes(amber.size(), eef1.size()));
		}

		@Override
		public void writeStaticStatic(MemoryBuffer buf) {

			ConfSpace.IndicesStatic amberIndices = confSpace.indicesStatic(0);
			ConfSpace.IndicesStatic eef1Indices = confSpace.indicesStatic(1);

			writeAtomPairs(
				new AtomPairWriter() {
					@Override public int size() { return amberIndices.size(); }
					@Override public int atomi1(int i) { return confSpace.getStaticAtomIndex(amberIndices.getStaticAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getStaticAtomIndex(amberIndices.getStaticAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(0, amberIndices.getParamsIndex(i)); }
				},
				new AtomPairWriter() {
					@Override public int size() { return eef1Indices.size(); }
					@Override public int atomi1(int i) { return confSpace.getStaticAtomIndex(eef1Indices.getStaticAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getStaticAtomIndex(eef1Indices.getStaticAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(1, eef1Indices.getParamsIndex(i)); }
				},
				buf
			);
		}

		@Override
		public void writeStaticPos(int posi1, int fragi1, MemoryBuffer buf) {

			ConfSpace.IndicesSingle amberIndices = confSpace.indicesSinglesByFrag(0, posi1, fragi1);
			ConfSpace.IndicesSingle eef1Indices = confSpace.indicesSinglesByFrag(1, posi1, fragi1);

			writeAtomPairs(
				new AtomPairWriter() {
					@Override public int size() { return amberIndices.sizeStatics(); }
					@Override public int atomi1(int i) { return confSpace.getStaticAtomIndex(amberIndices.getStaticStaticAtomIndex(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi1, amberIndices.getStaticConfAtomIndex(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(0, amberIndices.getStaticParamsIndex(i)); }
				},
				new AtomPairWriter() {
					@Override public int size() { return eef1Indices.sizeStatics(); }
					@Override public int atomi1(int i) { return confSpace.getStaticAtomIndex(eef1Indices.getStaticStaticAtomIndex(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi1, eef1Indices.getStaticConfAtomIndex(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(1, eef1Indices.getStaticParamsIndex(i)); }
				},
				buf
			);
		}

		@Override
		public void writePos(int posi1, int fragi1, MemoryBuffer buf) {

			ConfSpace.IndicesSingle amberIndices = confSpace.indicesSinglesByFrag(0, posi1, fragi1);
			ConfSpace.IndicesSingle eef1Indices = confSpace.indicesSinglesByFrag(1, posi1, fragi1);

			writeAtomPairs(
				new AtomPairWriter() {
					@Override public int size() { return amberIndices.sizeInternals(); }
					@Override public int atomi1(int i) { return confSpace.getConfAtomIndex(posi1, amberIndices.getInternalConfAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi1, amberIndices.getInternalConfAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(0, amberIndices.getInternalParamsIndex(i)); }
				},
				new AtomPairWriter() {
					@Override public int size() { return eef1Indices.sizeInternals(); }
					@Override public int atomi1(int i) { return confSpace.getConfAtomIndex(posi1, eef1Indices.getInternalConfAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi1, eef1Indices.getInternalConfAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(1, eef1Indices.getInternalParamsIndex(i)); }
				},
				buf
			);
		}

		@Override
		public void writePosPos(int posi1, int fragi1, int posi2, int fragi2, MemoryBuffer buf) {

			ConfSpace.IndicesPair amberIndices = confSpace.indicesPairsByFrags(0, posi1, fragi1, posi2, fragi2);
			ConfSpace.IndicesPair eef1Indices = confSpace.indicesPairsByFrags(1, posi1, fragi1, posi2, fragi2);

			writeAtomPairs(
				new AtomPairWriter() {
					@Override public int size() { return amberIndices.size(); }
					@Override public int atomi1(int i) { return confSpace.getConfAtomIndex(posi1, amberIndices.getConfAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi2, amberIndices.getConfAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(0, amberIndices.getParamsIndex(i)); }
				},
				new AtomPairWriter() {
					@Override public int size() { return eef1Indices.size(); }
					@Override public int atomi1(int i) { return confSpace.getConfAtomIndex(posi1, eef1Indices.getConfAtom1Index(i)); }
					@Override public int atomi2(int i) { return confSpace.getConfAtomIndex(posi2, eef1Indices.getConfAtom2Index(i)); }
					@Override public double[] params(int i) { return confSpace.ffparams(1, eef1Indices.getParamsIndex(i)); }
				},
				buf
			);
		}
	}

	public final ConfSpace confSpace;
	public final Precision precision;
	public final ForcefieldsImpl forcefieldsImpl;

	private final MemoryBuffer buf;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	class SConfSpace extends Struct {
		Int32 num_pos = int32();
		Int32 max_num_conf_atoms = int32();
		Int32 max_num_dofs = int32();
		Int32 num_molecule_motions = int32();
		Int64 positions_offset = int64();
		Int64 static_atom_coords_offset = int64();
		Int64 static_atom_molis_offset = int64();
		Int64 params_offset = int64();
		Int64 pos_pairs_offset = int64();
		Int64 molecule_motions_offset = int64();
		Real static_energy;
		Pad pad;
		void init() {
			static_energy = real(precision);
			pad = pad(precision.map(4, 0));
			init(
				precision.map(72, 72),
				"num_pos", "max_num_conf_atoms", "max_num_dofs", "num_molecule_motions",
				"positions_offset", "static_atom_coords_offset", "static_atom_molis_offset", "params_offset", "pos_pairs_offset",
				"molecule_motions_offset",
				"static_energy", "pad"
			);
		}
	}
	private final SConfSpace confSpaceStruct = new SConfSpace();

	static class SPos extends Struct {
		Int32 num_confs = int32();
		Int32 max_num_atoms = int32();
		Int32 num_frags = int32();
		Pad pad;
		void init() {
			pad = pad(4);
			init(
				16, "num_confs", "max_num_atoms", "num_frags", "pad"
			);
		}
	}
	private final SPos posStruct = new SPos();

	class SConf extends Struct {
		Int64 atom_coords_offset = int64();
		Int64 atom_molis_offset = int64();
		Int32 frag_index = int32();
		Pad pad;
		Real internal_energy;
		Int64 num_motions = int64();
		Int64 motions_offset = int64();
		void init() {
			pad = pad(precision.map(0, 4));
			internal_energy = real(precision);
			init(
				precision.map(40, 48),
				"atom_coords_offset", "atom_molis_offset", "frag_index", "pad", "internal_energy",
				"num_motions", "motions_offset"
			);
		}
	}
	private final SConf confStruct = new SConf();

	static class SArray extends Struct {
		Int64 size = int64();
		Int64 things_ptr = int64();
		void init() {
			init(
				16, "size", "things_ptr"
			);
		}
	}
	private final SArray arrayStruct = new SArray();

	class SReal3 extends Struct {
		Real x;
		Real y;
		Real z;
		Pad pad;
		void init() {
			x = real(precision);
			y = real(precision);
			z = real(precision);
			pad = pad(precision.map(4, 0));
			init(
				precision.map(16, 24),
				"x", "y", "z", "pad"
			);
		}
	}
	private final SReal3 real3Struct = new SReal3();

	class SPosInter extends Struct {
		Int32 posi1 = int32();
		Int32 posi2 = int32();
		Real weight;
		Real offset;
		void init() {
			weight = real(precision);
			offset = real(precision);
			init(
				precision.map(16, 24),
				"posi1", "posi2", "weight", "offset"
			);
		}
	}
	private final SPosInter posInterStruct = new SPosInter();

	class SDihedral extends Struct {
		Real min_radians;
		Real max_radians;
		Int32 a_index = int32();
		Int32 b_index = int32();
		Int32 c_index = int32();
		Int32 d_index = int32();
		Int32 num_rotated = int32();
		Int32 modified_posi = int32();

		void init() {
			min_radians = real(precision);
			max_radians = real(precision);
			init(
				precision.map(32, 40),
				"min_radians", "max_radians",
				"a_index", "b_index", "c_index", "d_index", "num_rotated",
				"modified_posi"
			);
		}

		long bytes(int numRotated) {
			return bytes() + MemoryBuffer.padToAlignment(Int32.bytes*numRotated, 8);
		}
	}
	private final SDihedral dihedralStruct = new SDihedral();
	private static final int dihedralId = 0;

	class STranslationRotation extends Struct {
		Real max_distance;
		Real max_radians;
		StructField<SReal3> centroid;
		Int32 moli = int32();
		Pad pad = pad(4);
		void init() {
			max_distance = real(precision);
			max_radians = real(precision);
			centroid = struct(real3Struct);
			init(
				precision.map(32, 48),
				"max_distance", "max_radians", "centroid",
				"moli", "pad"
			);
		}
	}
	private final STranslationRotation transRotStruct = new STranslationRotation();
	private static final int transRotId = 1;

	public NativeConfEnergyCalculator(ConfSpace confSpace, Precision precision) {

		this.confSpace = confSpace;
		this.precision = precision;

		// find the forcefield implementation, or die trying
		EnergyCalculator.Type[] ecalcTypes = Arrays.stream(confSpace.ecalcs)
			.map(EnergyCalculator::type)
			.toArray(EnergyCalculator.Type[]::new);
		if (ecalcTypes.length == 2 && ecalcTypes[0] == AmberEnergyCalculator.type && ecalcTypes[1] == EEF1EnergyCalculator.type) {
			forcefieldsImpl = new AmberEef1();
		} else {
			throw new IllegalArgumentException("No native implementation for forcefields: " + Arrays.toString(ecalcTypes));
		}

		// once we know the precision, init the structs
		confSpaceStruct.init();
		posStruct.init();
		confStruct.init();
		arrayStruct.init();
		real3Struct.init();
		posInterStruct.init();
		dihedralStruct.init();
		transRotStruct.init();

		Int64.Array posOffsets = int64array();
		Int64.Array confOffsets = int64array();
		Int64.Array posPairOffsets = int64array();
		Int64.Array fragOffsets = int64array();
		Int64.Array molMotionOffsets = int64array();
		Int64.Array confMotionOffsets = int64array();

		long motionIdBytes = Int64.bytes;

		// calculate how much memory we need for the conf space buffer
		long bufSize = confSpaceStruct.bytes()
			+ posOffsets.bytes(confSpace.positions.length)
			+ sum(confSpace.positions, pos ->
				posStruct.bytes()
				+ confOffsets.bytes(pos.confs.length)
				+ sum(pos.confs, conf ->
					confStruct.bytes()
					+ arrayStruct.bytes()
					+ real3Struct.bytes()*conf.coords.size
					+ arrayStruct.bytes()
					+ MemoryBuffer.padToAlignment(Int32.bytes*conf.coords.size, 8)
					+ confMotionOffsets.bytes(conf.motions.length)
					+ motionIdBytes*conf.motions.length
					+ sum(conf.motions, motion -> {
						if (motion instanceof DihedralAngle.Description) {
							var dihedral = (DihedralAngle.Description)motion;
							return dihedralStruct.bytes(dihedral.rotated.length);
						} else {
							throw new UnsupportedOperationException(motion.getClass().getName());
						}
					})
				)
			)
			+ arrayStruct.bytes()
			+ real3Struct.bytes()*confSpace.staticCoords.size
			+ arrayStruct.bytes()
			+ MemoryBuffer.padToAlignment(Int32.bytes*confSpace.staticCoords.size, 8)
			+ forcefieldsImpl.paramsBytes();

		// add space for the pos pair offsets
		int numPosPairs =
			1 // static-static
			+ confSpace.numPos() // static-pos
			+ confSpace.numPos() // pos
			+ confSpace.numPos()*(confSpace.numPos() - 1)/2; // pos-pos
		bufSize += posPairOffsets.bytes(numPosPairs);

		// add space for the static-static pairs
		bufSize += forcefieldsImpl.staticStaticBytes();

		// add space for the static-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.staticPosBytes(posi1, fragi1);
			}
		}

		// add space for the pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.posBytes(posi1, fragi1);
			}
		}

		// add space for the pos-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1)*confSpace.numFrag(posi2);
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						bufSize += forcefieldsImpl.posPosBytes( posi1, fragi1, posi2, fragi2);
					}
				}
			}
		}

		// add space for the molecule motions
		int numMolMotions = 0;
		for (int moli=0; moli<confSpace.molInfos.length; moli++) {
			var molInfo = confSpace.molInfos[moli];
			for (var motion : molInfo.motions) {
				numMolMotions += 1;
				if (motion instanceof DihedralAngle.Description) {
					var dihedral = (DihedralAngle.Description)motion;
					bufSize += dihedralStruct.bytes(dihedral.rotated.length);
				} else if (motion instanceof TranslationRotation.Description) {
					bufSize += transRotStruct.bytes();
				} else {
					throw new UnsupportedOperationException(motion.getClass().getName());
				}
			}
		}
		bufSize += molMotionOffsets.bytes(numMolMotions)
			+ motionIdBytes*numMolMotions;

		// allocate the buffer for the conf space
		buf = new MemoryBuffer(bufSize);

		// write the header
		var confSpaceAddr = buf.place(confSpaceStruct);
		confSpaceStruct.num_pos.set(confSpaceAddr, confSpace.positions.length);
		confSpaceStruct.max_num_conf_atoms.set(confSpaceAddr, confSpace.maxNumConfAtoms);
		confSpaceStruct.max_num_dofs.set(confSpaceAddr, confSpace.maxNumDofs);
		confSpaceStruct.num_molecule_motions.set(confSpaceAddr, numMolMotions);
		// we'll go back and write the offsets later
		confSpaceStruct.static_energy.set(confSpaceAddr, Arrays.stream(confSpace.staticEnergies).sum());

		// leave space for the position offsets
		confSpaceStruct.positions_offset.set(confSpaceAddr, buf.getPos());
		var posOffsetsAddr = buf.place(posOffsets, confSpace.positions.length);

		// write the positions
		for (ConfSpace.Pos pos : confSpace.positions) {
			posOffsets.set(posOffsetsAddr, pos.index, buf.getPos());
			var posAddr = buf.place(posStruct);
			posStruct.num_confs.set(posAddr, pos.confs.length);
			posStruct.max_num_atoms.set(posAddr, pos.maxNumAtoms);
			posStruct.num_frags.set(posAddr, pos.numFrags);

			// put the conf offsets
			var confOffsetsAddr = buf.place(confOffsets, pos.confs.length);
			// we'll go back and write them later though

			// write the confs
			for (ConfSpace.Conf conf : pos.confs) {

				confOffsets.set(confOffsetsAddr, conf.index, buf.getPos());
				var confAddr = buf.place(confStruct);
				confStruct.frag_index.set(confAddr, conf.fragIndex);
				confStruct.internal_energy.set(confAddr, Arrays.stream(conf.energies).sum());
				confStruct.num_motions.set(confAddr, conf.motions.length);

				// write the atom coords
				confStruct.atom_coords_offset.set(confAddr, buf.getPos());
				var atomCoordsAddr = buf.place(arrayStruct);
				arrayStruct.size.set(atomCoordsAddr, conf.coords.size);
				arrayStruct.things_ptr.set(atomCoordsAddr, 0);
				for (int i=0; i<conf.coords.size; i++) {
					var addr = buf.place(real3Struct);
					real3Struct.x.set(addr, conf.coords.x(i));
					real3Struct.y.set(addr, conf.coords.y(i));
					real3Struct.z.set(addr, conf.coords.z(i));
				}

				// write the atom molecule indices
				confStruct.atom_molis_offset.set(confAddr, buf.getPos());
				var atomMolisAddr = buf.place(arrayStruct);
				arrayStruct.size.set(atomMolisAddr, conf.coords.size);
				arrayStruct.things_ptr.set(atomMolisAddr, 0);
				for (int i=0; i<conf.numAtoms; i++) {
					buf.int32(conf.atomMolInfoIndices[i]);
				}
				buf.skipToAlignment(8);

				// write the motions
				confStruct.motions_offset.set(confAddr, buf.getPos());
				var confMotionAddr = buf.place(confMotionOffsets, conf.motions.length);
				for (int i=0; i<conf.motions.length; i++) {
					var motion = conf.motions[i];
					confMotionOffsets.set(confMotionAddr, i, buf.getPos());
					if (motion instanceof DihedralAngle.Description) {
						writeDihedral((DihedralAngle.Description)motion, pos.index, buf);
					} else {
						throw new UnsupportedOperationException(motion.getClass().getName());
					}
				}
			}
		}

		// write the static atom coords
		confSpaceStruct.static_atom_coords_offset.set(confSpaceAddr, buf.getPos());
		var atomCoordsAddr = buf.place(arrayStruct);
		arrayStruct.size.set(atomCoordsAddr, confSpace.staticCoords.size);
		arrayStruct.things_ptr.set(atomCoordsAddr, 0);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			var addr = buf.place(real3Struct);
			real3Struct.x.set(addr, confSpace.staticCoords.x(i));
			real3Struct.y.set(addr, confSpace.staticCoords.y(i));
			real3Struct.z.set(addr, confSpace.staticCoords.z(i));
		}

		// write the static atom molecule indices
		confSpaceStruct.static_atom_molis_offset.set(confSpaceAddr, buf.getPos());
		var atomMolisAddr = buf.place(arrayStruct);
		arrayStruct.size.set(atomMolisAddr, confSpace.staticCoords.size);
		arrayStruct.things_ptr.set(atomMolisAddr, 0);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			buf.int32(confSpace.staticMolInfoIndices[i]);
		}
		buf.skipToAlignment(8);

		// write the forcefield params
		confSpaceStruct.params_offset.set(confSpaceAddr, buf.getPos());
		forcefieldsImpl.writeParams(buf);

		// write the pos pairs
		confSpaceStruct.pos_pairs_offset.set(confSpaceAddr, buf.getPos());
		var posPairOffsetsAddr = buf.place(posPairOffsets, numPosPairs);

		// write the static-static pair
		posPairOffsets.set(posPairOffsetsAddr, 0, buf.getPos());
		forcefieldsImpl.writeStaticStatic(buf);

		// write the static-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(posPairOffsetsAddr, 1 + posi1, buf.getPos());
			var fragOffsetsAddr = buf.place(fragOffsets, confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragOffsetsAddr, fragi1, buf.getPos());
				forcefieldsImpl.writeStaticPos(posi1, fragi1, buf);
			}
		}

		// write the pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(posPairOffsetsAddr, 1 + confSpace.positions.length + posi1, buf.getPos());
			var fragOffsetsAddr = buf.place(fragOffsets, confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragOffsetsAddr, fragi1, buf.getPos());
				forcefieldsImpl.writePos(posi1, fragi1, buf);
			}
		}

		// write the pos-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				posPairOffsets.set(posPairOffsetsAddr, 1 + 2*confSpace.positions.length + posi1*(posi1 - 1)/2 + posi2, buf.getPos());
				var fragOffsetsAddr = buf.place(fragOffsets, confSpace.numFrag(posi1)*confSpace.numFrag(posi2));
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						fragOffsets.set(fragOffsetsAddr, fragi1*confSpace.numFrag(posi2) + fragi2, buf.getPos());
						forcefieldsImpl.writePosPos(posi1, fragi1, posi2, fragi2, buf);
					}
				}
			}
		}

		// write the molecule motions
		confSpaceStruct.molecule_motions_offset.set(confSpaceAddr, buf.getPos());
		var molMotionOffsetsAddr = buf.place(molMotionOffsets, numMolMotions);
		int molMotionIndex = 0;
		for (int moli=0; moli<confSpace.molInfos.length; moli++) {
			var molInfo = confSpace.molInfos[moli];
			for (var motion : molInfo.motions) {
				molMotionOffsets.set(molMotionOffsetsAddr, molMotionIndex++, buf.getPos());
				if (motion instanceof DihedralAngle.Description) {
					writeDihedral((DihedralAngle.Description)motion, PosInter.StaticPos, buf);
				} else if (motion instanceof TranslationRotation.Description) {
					writeTranslationRotation((TranslationRotation.Description)motion, moli, buf);
				} else {
					throw new UnsupportedOperationException(motion.getClass().getName());
				}
			}
		}

		assert(buf.getPos() == bufSize) : String.format("%d bytes leftover", bufSize - buf.getPos());
	}

	private void writeDihedral(DihedralAngle.Description desc, int posi, MemoryBuffer buf) {

		// write the motion id
		buf.int64(dihedralId);

		var dihedralAddr = buf.place(dihedralStruct);
		dihedralStruct.min_radians.set(dihedralAddr, Math.toRadians(desc.minDegrees));
		dihedralStruct.max_radians.set(dihedralAddr, Math.toRadians(desc.maxDegrees));
		dihedralStruct.a_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.a));
		dihedralStruct.b_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.b));
		dihedralStruct.c_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.c));
		dihedralStruct.d_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.d));
		dihedralStruct.num_rotated.set(dihedralAddr, desc.rotated.length);
		dihedralStruct.modified_posi.set(dihedralAddr, posi);

		var rotatedIndices = int32array();
		var indicesAddr = buf.place(rotatedIndices, desc.rotated.length);
		for (int i=0; i<desc.rotated.length; i++) {
			rotatedIndices.set(indicesAddr, i, desc.getAtomIndex(confSpace, posi, desc.rotated[i]));
		}

		buf.skipToAlignment(8);
	}

	private void writeTranslationRotation(TranslationRotation.Description desc, int moli, MemoryBuffer buf) {

		// write the motion id
		buf.int64(transRotId);

		var transrotAddr = buf.place(transRotStruct);
		transRotStruct.max_distance.set(transrotAddr, desc.maxDistance);
		transRotStruct.max_radians.set(transrotAddr, desc.maxRotationRadians);
		var centroidAddr = transRotStruct.centroid.offsetOf(transrotAddr);
		real3Struct.x.set(centroidAddr, desc.centroid.x);
		real3Struct.y.set(centroidAddr, desc.centroid.y);
		real3Struct.z.set(centroidAddr, desc.centroid.z);
		transRotStruct.moli.set(transrotAddr, moli);
	}

	@Override
	public Precision precision() {
		return precision;
	}

	public String version() {
		return String.format("%d.%d", NativeLib.version_major(), NativeLib.version_minor());
	}

	public AssignedCoords assign(int[] conf) {
		try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
			var confSpaceMem = this.buf;
			switch (precision) {
				case Float32 -> NativeLib.assign_f32(confSpaceMem.asByteBuffer(), conf, coordsMem.asByteBuffer());
				case Float64 -> NativeLib.assign_f64(confSpaceMem.asByteBuffer(), conf, coordsMem.asByteBuffer());
			}
			return makeCoords(coordsMem, conf);
		}
	}

	private AssignedCoords makeCoords(MemoryBuffer mem, int[] assignments) {

		AssignedCoords coords = new AssignedCoords(confSpace, assignments);

		// copy the coords from the native memory
		var arrayAddr = getArrayAddress(mem);
		for (int i=0; i<confSpace.maxNumConfAtoms; i++) {
			var addr = arrayAddr.sliceFrom((i*real3Struct.bytes()));
			coords.coords.set(
				i,
				real3Struct.x.get(addr),
				real3Struct.y.get(addr),
				real3Struct.z.get(addr)
			);
		}

		return coords;
	}

	@Override
	public void close() {
		buf.close();
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {
		try (var intersMem = makeIntersMem(inters)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				var confSpaceMem = this.buf;
				double energy = forcefieldsImpl.calc(confSpaceMem.asByteBuffer(), conf, intersMem.asByteBuffer(), coordsMem.asByteBuffer());
				return new EnergiedCoords(
					makeCoords(coordsMem, conf),
					energy
				);
			}
		}
	}

	@Override
	public double calcEnergy(int[] conf, List<PosInter> inters) {
		try (var intersMem = makeIntersMem(inters)) {
			var confSpaceMem = this.buf;
			return forcefieldsImpl.calc(confSpaceMem.asByteBuffer(), conf, intersMem.asByteBuffer(), null);
		}
	}

	private DoubleMatrix1D makeDofs(MemoryBuffer mem) {

		int size = (int)getArraySize(mem);
		DoubleMatrix1D vals = DoubleFactory1D.dense.make(size);

		Real.Array floats = realarray(precision);
		var addr = getArrayAddress(mem);
		for (int i=0; i<size; i++) {
			vals.set(i, floats.get(addr, i));
		}

		return vals;
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {
		try (var intersMem = makeIntersMem(inters)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				try (var dofsMem = makeArray(confSpace.maxNumDofs, precision.bytes)) {
					double energy;
					var confSpaceMem = this.buf;
					energy = forcefieldsImpl.minimize(
						confSpaceMem.asByteBuffer(), conf,
						intersMem.asByteBuffer(),
						coordsMem.asByteBuffer(), dofsMem.asByteBuffer()
					);
					return new EnergiedCoords(
						makeCoords(coordsMem, conf),
						energy,
						makeDofs(dofsMem)
					);
				}
			}
		}
	}

	@Override
	public double minimizeEnergy(int[] conf, List<PosInter> inters) {
		try (var intersMem = makeIntersMem(inters)) {
			var confSpaceMem = this.buf;
			return forcefieldsImpl.minimize(
				confSpaceMem.asByteBuffer(), conf,
				intersMem.asByteBuffer(),
				null, null
			);
		}
	}

	private MemoryBuffer makeIntersMem(List<PosInter> inters) {
		var buf = makeArray(inters.size(), posInterStruct.bytes());
		var arrayBuf = getArrayAddress(buf);
		for (var inter : inters) {
			var addr = arrayBuf.place(posInterStruct);
			posInterStruct.posi1.set(addr, inter.posi1);
			posInterStruct.posi2.set(addr, inter.posi2);
			posInterStruct.weight.set(addr, inter.weight);
			posInterStruct.offset.set(addr, inter.offset);
		}
		return buf;
	}

	// helpers for the Array class on the c++ size

	private MemoryBuffer makeArray(long size, long itemBytes) {
		MemoryBuffer buf = new MemoryBuffer(Int64.bytes*2 + size*itemBytes);
		buf.int64(size);
		buf.int64(0); // write 0 here so the c++ side knows the array came from the Java side
		return buf.sliceFrom(0); // this effectively rewinds the position counter
	}

	private long getArraySize(MemoryBuffer mem) {
		return mem.getLong(0);
	}

	private MemoryBuffer getArrayAddress(MemoryBuffer mem) {
		return mem.sliceFrom(Int64.bytes*2);
	}
}
