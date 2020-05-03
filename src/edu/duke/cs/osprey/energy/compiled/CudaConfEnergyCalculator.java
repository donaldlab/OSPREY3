package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.gpu.BufWriter;
import jdk.incubator.foreign.MemoryAddress;
import jdk.incubator.foreign.MemoryHandles;
import jdk.incubator.foreign.MemorySegment;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import static edu.duke.cs.osprey.gpu.Structs.*;
import static edu.duke.cs.osprey.tools.Log.log;


public class CudaConfEnergyCalculator implements ConfEnergyCalculator {

	private static class NativeLib {

		static {
			Native.register("CudaConfEcalc");
		}

		public static native int version_major();
		public static native int version_minor();
		public static native int cuda_version_driver();
		public static native int cuda_version_runtime();
		public static native int cuda_version_required();

		public static native Pointer alloc_conf_space_f32(ByteBuffer conf_space);
		public static native Pointer alloc_conf_space_f64(ByteBuffer conf_space);
		public static native void free_conf_space(Pointer p);

		public static native void assign_f32(Pointer conf_space, ByteBuffer conf, ByteBuffer out);
		public static native void assign_f64(Pointer conf_space, ByteBuffer conf, ByteBuffer out);
		public static native float calc_amber_eef1_f32(Pointer conf_space, ByteBuffer conf, ByteBuffer inters, ByteBuffer out_coords, long num_atoms);
		public static native double calc_amber_eef1_f64(Pointer conf_space, ByteBuffer conf, ByteBuffer inters, ByteBuffer out_coords, long num_atoms);
	}

	public static boolean isSupported() {
		return NativeLib.cuda_version_driver() > 0 && NativeLib.cuda_version_runtime() > 0;
	}

	/**
	 * Throws a helpful error message if this energy calculator is not supported.
	 */
	public static void checkSupported() {

		Function<Integer,String> versionString = v -> String.format("%d.%d", v/1000, (v % 1000)/10);

		int vDriver = NativeLib.cuda_version_driver();
		int vRequired = NativeLib.cuda_version_required();
		int vRuntime = NativeLib.cuda_version_runtime();

		if (vDriver <= 0) {
			throw new RuntimeException("No CUDA driver installed");
		}

		switch (vRuntime) {
			case -1: throw new RuntimeException("CUDA driver is insufficient."
				+ " Driver supports CUDA " + versionString.apply(vDriver)
				+ ", but CUDA " + versionString.apply(vRequired) + " is needed.");
			case -2: throw new RuntimeException("No CUDA device. Does this machine have an Nvidia GPU?");
			case Integer.MIN_VALUE: throw new RuntimeException("Unrecognized error: " + vRuntime);
			default: break;
		}
	}

	private interface ForcefieldsImpl {
		double calc(Pointer pConfSpace, ByteBuffer confBuf, ByteBuffer intersBuf, ByteBuffer coordsBuf, long numAtoms);
		double minimize(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords, ByteBuffer dofs);
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

			void setParams(double[] params) {
				esQ.set(params[0]);
				vdwA.set(params[1]);
				vdwB.set(params[2]);
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

			void setParams(double[] params) {
				vdwRadius1.set(params[0]);
				lambda1.set(params[1]);
				vdwRadius2.set(params[2]);
				lambda2.set(params[3]);
				alpha1.set(params[4]);
				alpha2.set(params[5]);
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
		public double calc(Pointer pConfSpace, ByteBuffer confBuf, ByteBuffer intersBuf, ByteBuffer coordsBuf, long numAtoms) {
			return switch (precision) {
				case Float32 -> NativeLib.calc_amber_eef1_f32(pConfSpace, confBuf, intersBuf, coordsBuf, numAtoms);
				case Float64 -> NativeLib.calc_amber_eef1_f64(pConfSpace, confBuf, intersBuf, coordsBuf, numAtoms);
			};
		}

		@Override
		public double minimize(ByteBuffer confSpaceBuf, int[] conf, ByteBuffer intersBuf, ByteBuffer coords, ByteBuffer dofs) {
			/*
			return switch (precision) {
				case Float32 -> NativeLib.minimize_amber_eef1_f32(confSpaceBuf, conf, intersBuf, coords, dofs);
				case Float64 -> NativeLib.minimize_amber_eef1_f64(confSpaceBuf, conf, intersBuf, coords, dofs);
			};
			*/
			throw new Error("TODO");
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

	private final Pointer pConfSpace;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	class SConfSpace extends Struct {
		Int32 num_pos = int32();
		Int32 max_num_conf_atoms = int32();
		Int32 max_num_dofs = int32();
		Int32 num_molecule_motions = int32();
		Int64 size = int64();
		Int64 positions_offset = int64();
		Int64 static_atoms_offset = int64();
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
				"size", "positions_offset", "static_atoms_offset", "params_offset", "pos_pairs_offset",
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
		Int64 atoms_offset = int64();
		Int32 frag_index = int32();
		Pad pad;
		Real internal_energy;
		Int64 num_motions = int64();
		Int64 motions_offset = int64();
		void init() {
			pad = pad(precision.map(0, 4));
			internal_energy = real(precision);
			init(
				precision.map(32, 40),
				"atoms_offset", "frag_index", "pad", "internal_energy",
				"num_motions", "motions_offset"
			);
		}
	}
	private final SConf confStruct = new SConf();

	static class SAtoms extends Struct {
		Int32 num_atoms = int32();
		Pad pad = pad(4);
		Int64 coords = int64();
		void init() {
			init(
				16, "num_atoms", "pad", "coords"
			);
		}
	}
	private final SAtoms atomsStruct = new SAtoms();

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
			return bytes() + BufWriter.padToAlignment(Int32.bytes*numRotated, 8);
		}
	}
	private final SDihedral dihedralStruct = new SDihedral();
	private static final int dihedralId = 0;

	class STranslationRotation extends Struct {
		Real max_distance;
		Real max_radians;
		StructField<SReal3> centroid;
		Int32 num_atoms = int32();
		Int32 num_modified_pos = int32();
		void init() {
			max_distance = real(precision);
			max_radians = real(precision);
			centroid = struct(real3Struct);
			init(
				precision.map(32, 48),
				"max_distance", "max_radians", "centroid",
				"num_atoms", "num_modified_pos"
			);
		}
	}
	private final STranslationRotation transRotStruct = new STranslationRotation();
	private static final int transRotId = 1;

	public CudaConfEnergyCalculator(ConfSpace confSpace, Precision precision) {

		this.confSpace = confSpace;
		this.precision = precision;

		checkSupported();

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
		atomsStruct.init();
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
					+ atomsStruct.bytes()
					+ real3Struct.bytes()*conf.coords.size
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
			bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.staticPosBytes(confSpace, posi1, fragi1);
			}
		}

		// add space for the pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				bufSize += forcefieldsImpl.posBytes(confSpace, posi1, fragi1);

			}
		}

		// add space for the pos-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				bufSize += fragOffsets.itemBytes*confSpace.numFrag(posi1)*confSpace.numFrag(posi2);
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						bufSize += forcefieldsImpl.posPosBytes(confSpace, posi1, fragi1, posi2, fragi2);
					}
				}
			}
		}

		// add space for the molecule motions
		int numMolMotions = 0;
		for (var molInfo : confSpace.molInfos) {
			for (var motion : molInfo.motions) {
				numMolMotions += 1;
				if (motion instanceof DihedralAngle.Description) {
					var dihedral = (DihedralAngle.Description)motion;
					bufSize += dihedralStruct.bytes(dihedral.rotated.length);
				} else {
					throw new UnsupportedOperationException(motion.getClass().getName());
				}
			}
		}
		bufSize += molMotionOffsets.bytes(numMolMotions)
			+ motionIdBytes*numMolMotions;

		// TODO: add space for the conf motions

		// allocate the buffer for the conf space
		try (MemorySegment confSpaceMem = MemorySegment.allocateNative(bufSize)) {
			BufWriter buf = new BufWriter(confSpaceMem);

			// write the header
			buf.place(confSpaceStruct);
			confSpaceStruct.num_pos.set(confSpace.positions.length);
			confSpaceStruct.max_num_conf_atoms.set(confSpace.maxNumConfAtoms);
			confSpaceStruct.max_num_dofs.set(confSpace.maxNumDofs);
			confSpaceStruct.num_molecule_motions.set(numMolMotions);
			confSpaceStruct.size.set(bufSize);
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
					confStruct.num_motions.set(conf.motions.length);

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

					// write the motions
					confStruct.motions_offset.set(buf.pos);
					buf.place(confMotionOffsets, conf.motions.length);
					for (int i=0; i<conf.motions.length; i++) {
						var motion = conf.motions[i];
						confMotionOffsets.set(i, buf.pos);
						if (motion instanceof DihedralAngle.Description) {
							writeDihedral((DihedralAngle.Description)motion, pos.index, buf);
						} else {
							throw new UnsupportedOperationException(motion.getClass().getName());
						}
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
				buf.place(fragOffsets, confSpace.numFrag(posi1));
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					fragOffsets.set(fragi1, buf.pos);
					forcefieldsImpl.writeStaticPos(confSpace, posi1, fragi1, buf);
				}
			}

			// write the pos pairs
			for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
				posPairOffsets.set(1 + confSpace.positions.length + posi1, buf.pos);
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
					buf.place(fragOffsets, confSpace.numFrag(posi1)*confSpace.numFrag(posi2));
					for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
						for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
							fragOffsets.set(fragi1*confSpace.numFrag(posi2) + fragi2, buf.pos);
							forcefieldsImpl.writePosPos(confSpace, posi1, fragi1, posi2, fragi2, buf);
						}
					}
				}
			}

			// write the molecule motions
			confSpaceStruct.molecule_motions_offset.set(buf.pos);
			buf.place(molMotionOffsets, numMolMotions);
			int molMotionIndex = 0;
			for (var molInfo : confSpace.molInfos) {
				for (var motion : molInfo.motions) {
					molMotionOffsets.set(molMotionIndex++, buf.pos);
					if (motion instanceof DihedralAngle.Description) {
						writeDihedral((DihedralAngle.Description)motion, PosInter.StaticPos, buf);
					} else {
						throw new UnsupportedOperationException(motion.getClass().getName());
					}
				}
			}

			assert(buf.pos == bufSize) : String.format("%d bytes leftover", bufSize - buf.pos);

			// upload the conf space to the GPU
			pConfSpace = switch (precision) {
				case Float32 -> NativeLib.alloc_conf_space_f32(confSpaceMem.asByteBuffer());
				case Float64 -> NativeLib.alloc_conf_space_f64(confSpaceMem.asByteBuffer());
			};
		}
	}

	private void writeDihedral(DihedralAngle.Description desc, int posi, BufWriter buf) {

		// write the motion id
		buf.int64(0);

		buf.place(dihedralStruct);
		dihedralStruct.min_radians.set(Math.toRadians(desc.minDegrees));
		dihedralStruct.max_radians.set(Math.toRadians(desc.maxDegrees));
		dihedralStruct.a_index.set(desc.getAtomIndex(confSpace, posi, desc.a));
		dihedralStruct.b_index.set(desc.getAtomIndex(confSpace, posi, desc.b));
		dihedralStruct.c_index.set(desc.getAtomIndex(confSpace, posi, desc.c));
		dihedralStruct.d_index.set(desc.getAtomIndex(confSpace, posi, desc.d));
		dihedralStruct.num_rotated.set(desc.rotated.length);
		dihedralStruct.modified_posi.set(posi);

		var rotatedIndices = int32array();
		buf.place(rotatedIndices, desc.rotated.length);
		for (int i=0; i<desc.rotated.length; i++) {
			rotatedIndices.set(i, desc.getAtomIndex(confSpace, posi, desc.rotated[i]));
		}

		buf.skipToAlignment(8);
	}

	@Override
	public Precision precision() {
		return precision;
	}

	public String version() {
		return String.format("%d.%d", NativeLib.version_major(), NativeLib.version_minor());
	}

	public AssignedCoords assign(int[] conf) {
		try (var confMem = makeConf(conf)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				switch (precision) {
					case Float32 -> NativeLib.assign_f32(pConfSpace, confMem.asByteBuffer(), coordsMem.asByteBuffer());
					case Float64 -> NativeLib.assign_f64(pConfSpace, confMem.asByteBuffer(), coordsMem.asByteBuffer());
				}
				return makeCoords(coordsMem, conf);
			}
		}
	}

	private MemorySegment makeConf(int[] conf) {

		MemorySegment mem = makeArray(confSpace.positions.length, Int32.bytes);
		var array = int32array();
		array.setAddress(getArrayAddress(mem));
		for (int posi=0; posi<confSpace.positions.length; posi++) {
			array.set(posi, conf[posi]);
		}

		return mem;
	}

	private AssignedCoords makeCoords(MemorySegment mem, int[] assignments) {

		AssignedCoords coords = new AssignedCoords(confSpace, assignments);

		// copy the coords from the native memory
		var addr = getArrayAddress(mem);
		for (int i=0; i<confSpace.maxNumConfAtoms; i++) {
			real3Struct.setAddress(addr.addOffset(i*real3Struct.bytes()));
			coords.coords.set(
				i,
				real3Struct.x.get(),
				real3Struct.y.get(),
				real3Struct.z.get()
			);
		}

		return coords;
	}

	@Override
	public void close() {
		NativeLib.free_conf_space(pConfSpace);
	}

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	@Override
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {
		try (var confMem = makeConf(conf)) {
			try (var intersMem = makeIntersMem(inters)) {
				try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
					double energy = forcefieldsImpl.calc(pConfSpace, confMem.asByteBuffer(), intersMem.asByteBuffer(), coordsMem.asByteBuffer(), confSpace.maxNumConfAtoms);
					return new EnergiedCoords(
						makeCoords(coordsMem, conf),
						energy
					);
				}
			}
		}
	}

	@Override
	public double calcEnergy(int[] conf, List<PosInter> inters) {
		try (var confMem = makeConf(conf)) {
			try (var intersMem = makeIntersMem(inters)) {
				return forcefieldsImpl.calc(pConfSpace, confMem.asByteBuffer(), intersMem.asByteBuffer(), null, confSpace.maxNumConfAtoms);
			}
		}
	}

	private DoubleMatrix1D makeDofs(MemorySegment mem) {

		int size = (int)getArraySize(mem);
		DoubleMatrix1D vals = DoubleFactory1D.dense.make(size);

		Real.Array floats = realarray(precision);
		floats.setAddress(getArrayAddress(mem));
		for (int i=0; i<size; i++) {
			vals.set(i, floats.get(i));
		}

		return vals;
	}

	@Override
	public EnergiedCoords minimize(int[] conf, List<PosInter> inters) {
		/* TODO
		try (var intersMem = makeIntersMem(inters)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				try (var dofsMem = makeArray(confSpace.maxNumDofs, precision.bytes)) {
					double energy;
					try (var confSpaceMem = this.confSpaceMem.acquire()) {
						energy = forcefieldsImpl.minimize(
							confSpaceMem.asByteBuffer(), conf,
							intersMem.asByteBuffer(),
							coordsMem.asByteBuffer(), dofsMem.asByteBuffer()
						);
					}
					return new EnergiedCoords(
						makeCoords(coordsMem, conf),
						energy,
						makeDofs(dofsMem)
					);
				}
			}
		}
		*/
		throw new Error("TODO");
	}

	@Override
	public double minimizeEnergy(int[] conf, List<PosInter> inters) {
		/* TODO
		try (var intersMem = makeIntersMem(inters)) {
			try (var confSpaceMem = this.confSpaceMem.acquire()) {
				return forcefieldsImpl.minimize(
					confSpaceMem.asByteBuffer(), conf,
					intersMem.asByteBuffer(),
					null, null
				);
			}
		}
		*/
		throw new Error("TODO");
	}

	private MemorySegment makeIntersMem(List<PosInter> inters) {
		MemorySegment mem = makeArray(inters.size(), posInterStruct.bytes());
		BufWriter buf = new BufWriter(mem);
		buf.pos = getArrayAddress(mem).offset();
		for (var inter : inters) {
			buf.place(posInterStruct);
			posInterStruct.posi1.set(inter.posi1);
			posInterStruct.posi2.set(inter.posi2);
			posInterStruct.weight.set(inter.weight);
			posInterStruct.offset.set(inter.offset);
		}
		return mem;
	}

	// helpers for the Array class on the c++ size

	private MemorySegment makeArray(long size, long itemBytes) {
		MemorySegment mem = MemorySegment.allocateNative(Int64.bytes*2 + size*itemBytes);
		BufWriter buf = new BufWriter(mem);
		buf.int64(size);
		buf.int64(0); // write 0 here so the c++ side knows the array came from the Java side
		return mem;
	}

	private long getArraySize(MemorySegment mem) {
		var h = MemoryHandles.varHandle(long.class, ByteOrder.nativeOrder());
		return (long)h.get(mem.baseAddress());
	}

	private MemoryAddress getArrayAddress(MemorySegment mem) {
		return mem.baseAddress().addOffset(Int64.bytes*2);
	}
}
