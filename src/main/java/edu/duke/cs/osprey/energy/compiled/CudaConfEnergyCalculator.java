package edu.duke.cs.osprey.energy.compiled;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import edu.duke.cs.osprey.confspace.compiled.AssignedCoords;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.motions.DihedralAngle;
import edu.duke.cs.osprey.confspace.compiled.motions.TranslationRotation;
import edu.duke.cs.osprey.gpu.MemoryBuffer;
import edu.duke.cs.osprey.gpu.Precision;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.gpu.Structs.*;
import static edu.duke.cs.osprey.tools.Log.log;


/**
 * An energy calculator that uses CUDA GPUs to accelerate the computation.
 */
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

		public static native Pointer alloc_gpu_infos();
		public static native void free_gpu_infos(Pointer p);

		public static native Pointer alloc_conf_space_f32(int device, ByteBuffer conf_space);
		public static native Pointer alloc_conf_space_f64(int device, ByteBuffer conf_space);
		public static native void free_conf_space(int device, Pointer p);

		public static native long minimize_batch_bufsize_device_f32(ByteBuffer conf_space_sizes, long batch_size);
		public static native long minimize_batch_bufsize_device_f64(ByteBuffer conf_space_sizes, long batch_size);
		public static native long minimize_batch_bufsize_host_f32(ByteBuffer conf_space_sizes, long batch_size);
		public static native long minimize_batch_bufsize_host_f64(ByteBuffer conf_space_sizes, long batch_size);

		public static native Pointer alloc_stream(int device, long host_bytes, long device_bytes);
		public static native void free_stream(int device, Pointer stream);

		public static native void assign_f32(int device, Pointer stream, Pointer conf_space, ByteBuffer conf, ByteBuffer out);
		public static native void assign_f64(int device, Pointer stream, Pointer conf_space, ByteBuffer conf, ByteBuffer out);
		public static native float calc_amber_eef1_f32(int device, Pointer stream, Pointer conf_space, ByteBuffer conf, ByteBuffer inters, ByteBuffer out_coords, long num_atoms);
		public static native double calc_amber_eef1_f64(int device, Pointer stream, Pointer conf_space, ByteBuffer conf, ByteBuffer inters, ByteBuffer out_coords, long num_atoms);
		public static native float minimize_amber_eef1_f32(int device, Pointer stream, Pointer conf_space, ByteBuffer conf_space_sizes, ByteBuffer jobs, ByteBuffer out_coords, ByteBuffer out_dofs);
		public static native double minimize_amber_eef1_f64(int device, Pointer stream, Pointer conf_space, ByteBuffer conf_space_sizes, ByteBuffer jobs, ByteBuffer out_coords, ByteBuffer out_dofs);
		public static native void minimize_batch_amber_eef1_f32(int device, Pointer stream, Pointer conf_space, ByteBuffer conf_space_sizes, ByteBuffer jobs, ByteBuffer out_energies);
		public static native void minimize_batch_amber_eef1_f64(int device, Pointer stream, Pointer conf_space, ByteBuffer conf_space_sizes, ByteBuffer jobs, ByteBuffer out_energies);
	}

	// TODO: could use records for this? It's a "preview" feature in JDK 14
	public static class GpuInfo {

		public final int id;
		public final String busId;
		public final String name;
		public final boolean integrated;
		public final boolean concurrentKernels;
		public final int numProcessors;
		public final int numAsyncEngines;
		public final long memTotal;
		public final long memFree;

		public GpuInfo(int id, String busId, String name, boolean integrated, boolean concurrentKernels, int numProcessors, int numAsyncEngines, long memTotal, long memFree) {
			this.id = id;
			this.busId = busId;
			this.name = name;
			this.integrated = integrated;
			this.concurrentKernels = concurrentKernels;
			this.numProcessors = numProcessors;
			this.numAsyncEngines = numAsyncEngines;
			this.memTotal = memTotal;
			this.memFree = memFree;
		}

		/** Returns the best number of streams/threads for saturating the device with work */
		public int bestNumStreams() {
			return 3;
		}

		/** Returns the best batch size for saturating the device with work */
		public int bestBatchSize() {
			return numProcessors/2;
		}

		@Override
		public int hashCode() {
			return busId.hashCode();
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof GpuInfo && equals((GpuInfo)other);
		}

		public boolean equals(GpuInfo other) {
			return this.busId.equals(other.busId);
		}

		@Override
		public String toString() {
			// this would be much easier with text block literals... sigh, at least it's a "preview" feature
			return "GPU[\n"
			+ String.format("%20s: %s\n", "id", id)
			+ String.format("%20s: %s\n", "bus id", busId)
			+ String.format("%20s: %s\n", "name", name)
			+ String.format("%20s: %b\n", "integrated", integrated)
			+ String.format("%20s: %b\n", "concurrent kernels", concurrentKernels)
			+ String.format("%20s: %d\n", "num processors", numProcessors)
			+ String.format("%20s: %d\n", "async engines", numAsyncEngines)
			+ String.format("%20s: %d B, %.1f GiB\n", "memory total", memTotal, memTotal/1024.0/1024.0/1024.0)
			+ String.format("%20s: %d B, %.1f GiB\n", "memory free", memFree, memFree/1024.0/1024.0/1024.0)
			+ "]";
		}
	}

	static class SArray extends Struct {
		Int64 size = int64();
		Int64 things_ptr = int64();
		SArray init() {
			init(16, "size", "things_ptr");
			return this;
		}
		long bytes(long numItems, long itemBytes) {
			return bytes() + padToGpuAlignment(numItems*itemBytes);
		}
	}
	private static final SArray arrayStruct = new SArray().init();

	static class SGpuInfo extends Struct {
		Pad bus_id = pad(16);
		Pad name = pad(256);
		Int32 integrated = int32();
		Int32 concurrent_kernels = int32();
		Int32 num_processors = int32();
		Int32 num_async_engines = int32();
		Int64 mem_total = int64();
		Int64 mem_free = int64();
		SGpuInfo init() {
			init(304,
				"bus_id", "name", "integrated", "concurrent_kernels",
				"num_processors", "num_async_engines", "mem_total", "mem_free"
			);
			return this;
		}
	}
	static final SGpuInfo gpuInfoStruct = new SGpuInfo().init();

	public static List<GpuInfo> getGpusInfos() {

		if (!isSupported()) {
			return Collections.emptyList();
		}

		Pointer pArray = NativeLib.alloc_gpu_infos();

		try {
			// convert the pointer to a byte buffer we can read in java-land
			int size = (int)pArray.getLong(0);
			ByteBuffer buf = pArray.getByteBuffer(0, arrayStruct.bytes(size, gpuInfoStruct.bytes()));

			List<GpuInfo> infos = new ArrayList<>(size);

			try (var mem = MemoryBuffer.ofByteBuffer(buf)) {

				var p = getArrayAddress(mem);
				var charsBusId = char8array();
				var charsName = char8array();

				for (int device=0; device<size; device++) {

					infos.add(new GpuInfo(
						device,
						charsBusId.getNullTerminated(p.sliceFrom(gpuInfoStruct.bus_id.offset()), (int)gpuInfoStruct.bus_id.bytes),
						charsName.getNullTerminated(p.sliceFrom(gpuInfoStruct.name.offset()), (int)gpuInfoStruct.name.bytes),
						gpuInfoStruct.integrated.get(p) != 0,
						gpuInfoStruct.concurrent_kernels.get(p) != 0,
						gpuInfoStruct.num_processors.get(p),
						gpuInfoStruct.num_async_engines.get(p),
						gpuInfoStruct.mem_total.get(p),
						gpuInfoStruct.mem_free.get(p)
					));

					p = p.sliceFrom(gpuInfoStruct.bytes());
				}

				return infos;
			}
		} finally {
			NativeLib.free_gpu_infos(pArray);
		}
	}

	// biggest LD instruction we can use, in bytes
	// to make sure structs get aligned so we can minimize LD instructions needed
	private static final long WidestGpuLoad = 16;

	private static long padToGpuAlignment(long pos) {
		return MemoryBuffer.padToAlignment(pos, WidestGpuLoad);
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
			case -3: throw new RuntimeException("CUDA failed to initialize");
			case -4: throw new RuntimeException(
				"This error indicates that the system was upgraded to run with forward compatibility"
				+ ", but the visible hardware detected by CUDA does not support this configuration."
				+ " Refer to the compatibility documentation for the supported hardware matrix or ensure"
				+ " that only supported hardware is visible during initialization via the CUDA_VISIBLE_DEVICES"
				+ " environment variable.");
			case Integer.MIN_VALUE: throw new RuntimeException("Unrecognized error");
			default: break;
		}
	}

	private interface ForcefieldsImpl {
		double calc(Stream stream, ByteBuffer confBuf, ByteBuffer intersBuf, ByteBuffer coordsBuf, long numAtoms);
		double minimize(Stream stream, ByteBuffer jobsBuf, ByteBuffer coordsBuf, ByteBuffer dofsBuf);
		void minimizeBatch(Stream stream, ByteBuffer jobsBuf, ByteBuffer energiesBuf);
		void dumpJobs(File file, ByteBuffer jobsBuf);
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
			Pad pad = pad(15);
			void init() {
				init(16, "distance_dependent_dielectric", "pad");
			}
		}
		final SParamsAmberEef1 paramsStruct = new SParamsAmberEef1();

		class SAtomPairs extends Struct {
			final Int32 num_amber = int32();
			final Int32 num_eef1 = int32();
			final Pad pad = pad(8);
			void init() {
				init(16, "num_amber", "num_eef1", "pad");
			}
		}
		final SAtomPairs atomPairsStruct = new SAtomPairs();

		class SAtomPairAmberF32a extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			void init() {
				init(8, "atomi1", "atomi2");
			}
		}
		final SAtomPairAmberF32a amberF32aStruct = new SAtomPairAmberF32a();

		class SAtomPairAmberF32b extends Struct {
			final Float32 esQ = float32();
			final Float32 vdwA = float32();
			final Float32 vdwB = float32();
			final Pad pad = pad(4);
			void init() {
				init(16, "esQ", "vdwA", "vdwB", "pad");
			}
			void setParams(MemoryBuffer buf, double[] params) {
				esQ.set(buf, (float)params[0]);
				vdwA.set(buf, (float)params[1]);
				vdwB.set(buf, (float)params[2]);
			}
		}
		final SAtomPairAmberF32b amberF32bStruct = new SAtomPairAmberF32b();

		class SAtomPairAmberF64a extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			final Float64 esQ = float64();
			void setParams(MemoryBuffer buf, double[] params) {
				esQ.set(buf, params[0]);
			}
			void init() {
				init(16, "atomi1", "atomi2", "esQ");
			}
		}
		final SAtomPairAmberF64a amberF64aStruct = new SAtomPairAmberF64a();

		class SAtomPairAmberF64b extends Struct {
			final Float64 vdwA = float64();
			final Float64 vdwB = float64();
			void init() {
				init(16, "vdwA", "vdwB");
			}
			void setParams(MemoryBuffer buf, double[] params) {
				vdwA.set(buf, params[1]);
				vdwB.set(buf, params[2]);
			}
		}
		final SAtomPairAmberF64b amberF64bStruct = new SAtomPairAmberF64b();

		class SAtomPairEef1F32a extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			void init() {
				init(8, "atomi1", "atomi2");
			}
		}
		final SAtomPairEef1F32a eef1F32aStruct = new SAtomPairEef1F32a();

		class SAtomPairEef1F32b extends Struct {
			final Float32 vdwRadius = float32();
			final Float32 oolambda = float32();
			final Float32 alpha = float32();
			final Pad pad = pad(4);
			void init() {
				init(16, "vdwRadius", "oolambda", "alpha", "pad");
			}
			void setParams1(MemoryBuffer buf, double[] params) {
				vdwRadius.set(buf, (float)params[0]);
				oolambda.set(buf, (float)(1.0/params[1]));
				alpha.set(buf, (float)params[4]);
			}
			void setParams2(MemoryBuffer buf, double[] params) {
				vdwRadius.set(buf, (float)params[2]);
				oolambda.set(buf, (float)(1.0/params[3]));
				alpha.set(buf, (float)params[5]);
			}
		}
		final SAtomPairEef1F32b eef1F32bStruct = new SAtomPairEef1F32b();

		class SAtomPairEef1F64a extends Struct {
			final Int32 atomi1 = int32();
			final Int32 atomi2 = int32();
			void init() {
				init(8, "atomi1", "atomi2");
			}
		}
		final SAtomPairEef1F64a eef1F64aStruct = new SAtomPairEef1F64a();

		class SAtomPairEef1F64b extends Struct {
			final Float64 vdwRadius = float64();
			final Float64 oolambda = float64();
			void init() {
				init(16, "vdwRadius", "oolambda");
			}
			void setParams1(MemoryBuffer buf, double[] params) {
				vdwRadius.set(buf, params[0]);
				oolambda.set(buf, 1.0/params[1]);
			}
			void setParams2(MemoryBuffer buf, double[] params) {
				vdwRadius.set(buf, params[2]);
				oolambda.set(buf, 1.0/params[3]);
			}
		}
		final SAtomPairEef1F64b eef1F64bStruct = new SAtomPairEef1F64b();

		class SAtomPairEef1F64c extends Struct {
			final Float64 alpha1 = float64();
			final Float64 alpha2 = float64();
			void init() {
				init(16, "alpha1", "alpha2");
			}
			void setParams(MemoryBuffer buf, double[] params) {
				alpha1.set(buf, params[4]);
				alpha2.set(buf, params[5]);
			}
		}
		final SAtomPairEef1F64c eef1F64cStruct = new SAtomPairEef1F64c();

		AmberEef1() {

			// once we know the precision, init the structs
			paramsStruct.init();
			atomPairsStruct.init();
			amberF32aStruct.init();
			amberF32bStruct.init();
			amberF64aStruct.init();
			amberF64bStruct.init();
			eef1F32aStruct.init();
			eef1F32bStruct.init();
			eef1F64aStruct.init();
			eef1F64bStruct.init();
			eef1F64cStruct.init();
		}

		@Override
		public double calc(Stream stream, ByteBuffer confBuf, ByteBuffer intersBuf, ByteBuffer coordsBuf, long numAtoms) {
			return switch (precision) {
				case Float32 -> (double)NativeLib.calc_amber_eef1_f32(stream.device, stream.stream, stream.pConfSpace, confBuf, intersBuf, coordsBuf, numAtoms);
				case Float64 -> NativeLib.calc_amber_eef1_f64(stream.device, stream.stream, stream.pConfSpace, confBuf, intersBuf, coordsBuf, numAtoms);
			};
		}

		@Override
		public double minimize(Stream stream, ByteBuffer jobsBuf, ByteBuffer coordsBuf, ByteBuffer dofsBuf) {
			return switch (precision) {
				case Float32 -> (double)NativeLib.minimize_amber_eef1_f32(stream.device, stream.stream, stream.pConfSpace, confSpaceSizesBuf, jobsBuf, coordsBuf, dofsBuf);
				case Float64 -> NativeLib.minimize_amber_eef1_f64(stream.device, stream.stream, stream.pConfSpace, confSpaceSizesBuf, jobsBuf, coordsBuf, dofsBuf);
			};
		}

		@Override
		public void minimizeBatch(Stream stream, ByteBuffer jobsBuf, ByteBuffer energiesBuf) {
			switch (precision) {
				case Float32 -> NativeLib.minimize_batch_amber_eef1_f32(stream.device, stream.stream, stream.pConfSpace, confSpaceSizesBuf, jobsBuf, energiesBuf);
				case Float64 -> NativeLib.minimize_batch_amber_eef1_f64(stream.device, stream.stream, stream.pConfSpace, confSpaceSizesBuf, jobsBuf, energiesBuf);
			}
		}

		@Override
		public void dumpJobs(File file, ByteBuffer jobsBuf) {

			try (var out = new DataOutputStream(new FileOutputStream(dumpJobsTo))) {

				log("dumping jobs ...");

				// write enough stuff to directly call minimize_batch_amber_eef1_fXX() on the C++ side
				out.writeInt(precision.ordinal());
				var spaceMem = confSpaceBuf;
				writeBuf(out, spaceMem.asByteBuffer());
				var sizesMem = confSpaceSizesMem;
				writeBuf(out, sizesMem.asByteBuffer());
				writeBuf(out, jobsBuf);

			} catch (IOException ex) {
				throw new Error("can't dump jobs", ex);
			}
		}

		private void writeBuf(DataOutputStream out, ByteBuffer buf)
		throws IOException {

			out.writeInt(buf.capacity());

			// this will likely be slow, but it's only for debugging, so performance here isn't super important
			for (int i=0; i<buf.capacity(); i++) {
				out.writeByte(buf.get(i));
			}
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
				+ switch (precision) {
					case Float32 ->
						padToGpuAlignment(numAmber*amberF32aStruct.bytes())
						+ numAmber*amberF32bStruct.bytes()
						+ padToGpuAlignment(numEef1*eef1F32aStruct.bytes())
						+ numEef1*eef1F32bStruct.bytes()
						+ numEef1*eef1F32bStruct.bytes();
					case Float64 ->
						padToGpuAlignment(numAmber*amberF64aStruct.bytes())
						+ numAmber*amberF64bStruct.bytes()
						+ padToGpuAlignment(numEef1*eef1F64aStruct.bytes())
						+ numEef1*eef1F64bStruct.bytes()
						+ numEef1*eef1F64bStruct.bytes()
						+ numEef1*eef1F64cStruct.bytes();
				};
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

			// make sure we're starting at optimal alignment for GPU loads
			assert (buf.isAligned(WidestGpuLoad)) : "Not aligned!";
			long firstPos = buf.getPos();

			var atomPairsAddr = buf.place(atomPairsStruct);
			atomPairsStruct.num_amber.set(atomPairsAddr, amber.size());
			atomPairsStruct.num_eef1.set(atomPairsAddr, eef1.size());

			switch (precision) {

				case Float32 -> {

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<amber.size(); i++) {
						var addr = buf.place(amberF32aStruct);
						int atomi1 = amber.atomi1(i);
						int atomi2 = amber.atomi2(i);
						assert (atomi1 != atomi2);
						amberF32aStruct.atomi1.set(addr, atomi1);
						amberF32aStruct.atomi2.set(addr, atomi2);
					}
					buf.skipToAlignment(WidestGpuLoad);

					for (int i=0; i<amber.size(); i++) {
						var addr = buf.place(amberF32bStruct);
						amberF32bStruct.setParams(addr, amber.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F32aStruct);
						int atomi1 = eef1.atomi1(i);
						int atomi2 = eef1.atomi2(i);
						assert (atomi1 != atomi2);
						eef1F32aStruct.atomi1.set(addr, atomi1);
						eef1F32aStruct.atomi2.set(addr, atomi2);
					}
					buf.skipToAlignment(WidestGpuLoad);

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F32bStruct);
						eef1F32bStruct.setParams1(addr, eef1.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F32bStruct);
						eef1F32bStruct.setParams2(addr, eef1.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));
				}

				case Float64 -> {

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<amber.size(); i++) {
						var addr = buf.place(amberF64aStruct);
						int atomi1 = amber.atomi1(i);
						int atomi2 = amber.atomi2(i);
						assert (atomi1 != atomi2);
						amberF64aStruct.atomi1.set(addr, atomi1);
						amberF64aStruct.atomi2.set(addr, atomi2);
						amberF64aStruct.setParams(addr, amber.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<amber.size(); i++) {
						var addr = buf.place(amberF64bStruct);
						amberF64bStruct.setParams(addr, amber.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F64aStruct);
						int atomi1 = eef1.atomi1(i);
						int atomi2 = eef1.atomi2(i);
						assert (atomi1 != atomi2);
						eef1F64aStruct.atomi1.set(addr, atomi1);
						eef1F64aStruct.atomi2.set(addr, atomi2);
					}
					buf.skipToAlignment(WidestGpuLoad);

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F64bStruct);
						eef1F64bStruct.setParams1(addr, eef1.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F64bStruct);
						eef1F64bStruct.setParams2(addr, eef1.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));

					for (int i=0; i<eef1.size(); i++) {
						var addr = buf.place(eef1F64cStruct);
						eef1F64cStruct.setParams(addr, eef1.params(i));
					}

					assert (buf.isAligned(WidestGpuLoad));
				}
			}


			// make sure we ended at optimal alignment for GPU loads
			assert (buf.getPos() - firstPos == atomPairsBytes(amber.size(), eef1.size()))
				: String.format("overshot by %d bytes", buf.getPos() - firstPos - atomPairsBytes(amber.size(), eef1.size()));
			assert (buf.isAligned(WidestGpuLoad)) : "Not aligned!";
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
	public final List<GpuStreams> gpuStreams;
	public final long maxBatchSize;
	public final ForcefieldsImpl forcefieldsImpl;

	/**
	 * Write minimization jobs to the file and throw an error instead of actually minimizing them.
	 * Useful for debugging GPU kernels without Java's extra complication.
	 */
	public CudaConfEnergyCalculator dumpJobsTo(File file) {
		dumpJobsTo = file;
		return this;
	}
	private File dumpJobsTo = null;

	private final MemoryBuffer confSpaceBuf;
	private final Map<GpuInfo,Pointer> pConfSpace = new HashMap<>();
	private final MemoryBuffer confSpaceSizesMem;
	private final ByteBuffer confSpaceSizesBuf;

	// NOTE: prefix the struct classes with S to avoid name collisions with the related Java classes

	class SConfSpace extends Struct {
		Int32 num_pos = int32();
		Int32 max_num_conf_atoms = int32();
		Int32 max_num_dofs = int32();
		Int32 num_molecule_motions = int32();
		Int64 size = int64();
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
				precision.map(80, 80),
				"num_pos", "max_num_conf_atoms", "max_num_dofs", "num_molecule_motions",
				"size", "positions_offset", "static_atom_coords_offset", "static_atom_molis_offset", "params_offset", "pos_pairs_offset",
				"molecule_motions_offset",
				"static_energy", "pad"
			);
		}
	}
	private final SConfSpace confSpaceStruct = new SConfSpace();

	static class SConfSpaceSizes extends Struct {
		Int64 num_pos = int64();
		Int64 max_num_inters = int64();
		Int64 num_atoms = int64();
		Int64 max_num_dofs = int64();
		Int64 num_mol_motions = int64();
		SConfSpaceSizes init() {
			init(40, "num_pos", "max_num_inters", "num_atoms", "max_num_dofs", "num_mol_motions");
			return this;
		}
	}
	private static final SConfSpaceSizes confSpaceSizesStruct = new SConfSpaceSizes().init();

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
		Pad pad1;
		Real internal_energy;
		Int64 num_motions = int64();
		Int64 motions_offset = int64();
		Pad pad2;
		void init() {
			pad1 = pad(precision.map(0, 4));
			internal_energy = real(precision);
			pad2 = pad(precision.map(8, 0));
			init(
				precision.map(48, 48),
				"atom_coords_offset", "atom_molis_offset", "frag_index", "pad1", "internal_energy",
				"num_motions", "motions_offset", "pad2"
			);
		}
	}
	private final SConf confStruct = new SConf();

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
		Int64 id = int64();
		Real min_radians;
		Real max_radians;
		Int32 a_index = int32();
		Int32 b_index = int32();
		Int32 c_index = int32();
		Int32 d_index = int32();
		Int32 num_rotated = int32();
		Int32 modified_posi = int32();
		Pad pad;

		void init() {
			min_radians = real(precision);
			max_radians = real(precision);
			pad = pad(precision.map(8, 0));
			init(
				precision.map(48, 48),
				"id", "min_radians", "max_radians",
				"a_index", "b_index", "c_index", "d_index", "num_rotated",
				"modified_posi", "pad"
			);
		}

		long bytes(int numRotated) {
			return bytes() + padToGpuAlignment(Int32.bytes*numRotated);
		}
	}
	private final SDihedral dihedralStruct = new SDihedral();
	private static final int dihedralId = 0;

	class STranslationRotation extends Struct {
		Int64 id = int64();
		Real max_distance;
		Real max_radians;
		StructField<SReal3> centroid;
		Int32 moli = int32();
		Pad pad = pad(12);
		void init() {
			max_distance = real(precision);
			max_radians = real(precision);
			centroid = struct(real3Struct);
			init(
				precision.map(48, 64),
				"id", "max_distance", "max_radians", "centroid",
				"moli", "pad"
			);
		}
	}
	private final STranslationRotation transRotStruct = new STranslationRotation();
	private static final int transRotId = 1;

	public static class GpuStreams {

		public final GpuInfo gpuInfo;
		public final int numStreams;

		public GpuStreams(GpuInfo gpuInfo, int numStreams) {
			this.gpuInfo = gpuInfo;
			this.numStreams = numStreams;
		}
	}

	/** use all the GPUs by default */
	public CudaConfEnergyCalculator(ConfSpace confSpace, Precision precision) {
		this(confSpace, precision, getGpusInfos());
	}

	/** get GPU info from the Parallelism instance */
	public CudaConfEnergyCalculator(ConfSpace confSpace, Precision precision, Parallelism parallelism) {
		this(confSpace, precision, getGpus(parallelism.numGpus));
	}

	private static List<GpuInfo> getGpus(int num) {

		var infos = getGpusInfos();

		if (num > infos.size()) {
			log("WARN: tried to pick %d GPU(s), but only %d GPU(s) were available", num, infos.size());
			num = infos.size();
		}

		// TODO: pick the "best" GPUs by some criterion?

		return infos.subList(0, num);
	}

	/** use just the given GPUs */
	public CudaConfEnergyCalculator(ConfSpace confSpace, Precision precision, List<GpuInfo> gpus) {
		this(
			confSpace,
			precision,
			gpus.stream()
				.map(gpu -> new CudaConfEnergyCalculator.GpuStreams(gpu, gpu.bestNumStreams()))
				.collect(Collectors.toList()),
			gpus.stream()
				.mapToInt(gpu -> gpu.bestBatchSize())
				.max()
				.orElse(Integer.MAX_VALUE)
		);
	}

	public CudaConfEnergyCalculator(ConfSpace confSpace, Precision precision, List<GpuStreams> gpuStreams, long maxBatchSize) {

		checkSupported();

		if (gpuStreams == null || gpuStreams.isEmpty()) {
			throw new IllegalArgumentException("0 GPUs selected");
		}

		this.confSpace = confSpace;
		this.precision = precision;
		this.gpuStreams = gpuStreams;
		this.maxBatchSize = maxBatchSize;

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
		long confSpaceSize = confSpaceStruct.bytes()
			+ padToGpuAlignment(posOffsets.bytes(confSpace.positions.length));

		long positionsSize = sum(confSpace.positions, pos ->
			posStruct.bytes()
			+ padToGpuAlignment(confOffsets.bytes(pos.confs.length))
			+ sum(pos.confs, conf ->
				confStruct.bytes()
				+ arrayStruct.bytes(conf.coords.size, real3Struct.bytes())
				+ arrayStruct.bytes(conf.coords.size, Int32.bytes)
				+ padToGpuAlignment(confMotionOffsets.bytes(conf.motions.length))
				+ sum(conf.motions, motion -> {
					if (motion instanceof DihedralAngle.Description) {
						var dihedral = (DihedralAngle.Description)motion;
						return dihedralStruct.bytes(dihedral.rotated.length);
					} else {
						throw new UnsupportedOperationException(motion.getClass().getName());
					}
				})
			)
		);

		long staticCoordsSize = arrayStruct.bytes(confSpace.staticCoords.size, real3Struct.bytes())
			+ arrayStruct.bytes(confSpace.staticCoords.size, Int32.bytes);

		long forcefieldSize = forcefieldsImpl.paramsBytes();

		// add space for the pos pair offsets
		int numPosPairs =
			1 // static-static
			+ confSpace.numPos() // static-pos
			+ confSpace.numPos() // pos
			+ confSpace.numPos()*(confSpace.numPos() - 1)/2; // pos-pos
		forcefieldSize += padToGpuAlignment(posPairOffsets.bytes(numPosPairs));

		// add space for the static-static pairs
		forcefieldSize += forcefieldsImpl.staticStaticBytes();

		// add space for the static-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			forcefieldSize += padToGpuAlignment(fragOffsets.itemBytes*confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				forcefieldSize += forcefieldsImpl.staticPosBytes(posi1, fragi1);
			}
		}

		// add space for the pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			forcefieldSize += padToGpuAlignment(fragOffsets.itemBytes*confSpace.numFrag(posi1));
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				forcefieldSize += forcefieldsImpl.posBytes(posi1, fragi1);

			}
		}

		// add space for the pos-pos pairs and offsets
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				forcefieldSize += padToGpuAlignment(fragOffsets.itemBytes*confSpace.numFrag(posi1)*confSpace.numFrag(posi2));
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						forcefieldSize += forcefieldsImpl.posPosBytes(posi1, fragi1, posi2, fragi2);
					}
				}
			}
		}

		// add space for molecule motions
		long molMotionsSize = 0;
		int numMolMotions = 0;
		for (int moli=0; moli<confSpace.molInfos.length; moli++) {
			var molInfo = confSpace.molInfos[moli];
			for (var motion : molInfo.motions) {
				numMolMotions += 1;
				if (motion instanceof DihedralAngle.Description) {
					var dihedral = (DihedralAngle.Description)motion;
					molMotionsSize += dihedralStruct.bytes(dihedral.rotated.length);
				} else if (motion instanceof TranslationRotation.Description) {
					molMotionsSize += transRotStruct.bytes();
				} else {
					throw new UnsupportedOperationException(motion.getClass().getName());
				}
			}
		}
		molMotionsSize += molMotionOffsets.bytes(numMolMotions);

		// allocate the buffer for the conf space
		long bufSize = confSpaceSize + positionsSize + staticCoordsSize + forcefieldSize + molMotionsSize;
		confSpaceBuf = new MemoryBuffer(bufSize);

		assert (confSpaceBuf.isAligned(WidestGpuLoad));

		// write the header
		var confSpaceAddr = confSpaceBuf.place(confSpaceStruct);
		confSpaceStruct.num_pos.set(confSpaceAddr, confSpace.positions.length);
		confSpaceStruct.max_num_conf_atoms.set(confSpaceAddr, confSpace.maxNumConfAtoms);
		confSpaceStruct.max_num_dofs.set(confSpaceAddr, confSpace.maxNumDofs);
		confSpaceStruct.num_molecule_motions.set(confSpaceAddr, numMolMotions);
		confSpaceStruct.size.set(confSpaceAddr, bufSize);
		// we'll go back and write the offsets later
		confSpaceStruct.static_energy.set(confSpaceAddr, Arrays.stream(confSpace.staticEnergies).sum());
		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		// leave space for the position offsets
		confSpaceStruct.positions_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		var posOffsetsAddr = confSpaceBuf.place(posOffsets, confSpace.positions.length, WidestGpuLoad);

		assert (confSpaceBuf.getPos() == confSpaceSize);
		assert (confSpaceBuf.isAligned(WidestGpuLoad));

		// write the positions
		for (ConfSpace.Pos pos : confSpace.positions) {

			assert (confSpaceBuf.isAligned(WidestGpuLoad));

			posOffsets.set(posOffsetsAddr, pos.index, confSpaceBuf.getPos());
			var posAddr = confSpaceBuf.place(posStruct);
			posStruct.num_confs.set(posAddr, pos.confs.length);
			posStruct.max_num_atoms.set(posAddr, pos.maxNumAtoms);
			posStruct.num_frags.set(posAddr, pos.numFrags);

			assert (confSpaceBuf.isAligned(WidestGpuLoad));

			// put the conf offsets
			var confOffsetsAddr = confSpaceBuf.place(confOffsets, pos.confs.length, WidestGpuLoad);
			// we'll go back and write them later though

			assert (confSpaceBuf.isAligned(WidestGpuLoad));

			// write the confs
			for (ConfSpace.Conf conf : pos.confs) {

				assert (confSpaceBuf.isAligned(WidestGpuLoad));

				confOffsets.set(confOffsetsAddr, conf.index, confSpaceBuf.getPos());
				var confAddr = confSpaceBuf.place(confStruct);
				confStruct.frag_index.set(confAddr, conf.fragIndex);
				confStruct.internal_energy.set(confAddr, Arrays.stream(conf.energies).sum());
				confStruct.num_motions.set(confAddr, conf.motions.length);
				confSpaceBuf.skipToAlignment(WidestGpuLoad);

				// write the atom coords
				confStruct.atom_coords_offset.set(confAddr, confSpaceBuf.getPos());
				var atomCoordsAddr = confSpaceBuf.place(arrayStruct);
				arrayStruct.size.set(atomCoordsAddr, conf.coords.size);
				arrayStruct.things_ptr.set(atomCoordsAddr, 0);
				for (int i=0; i<conf.coords.size; i++) {
					var addr = confSpaceBuf.place(real3Struct);
					real3Struct.x.set(addr, conf.coords.x(i));
					real3Struct.y.set(addr, conf.coords.y(i));
					real3Struct.z.set(addr, conf.coords.z(i));
				}
				confSpaceBuf.skipToAlignment(WidestGpuLoad);

				// write the atom molis
				confStruct.atom_molis_offset.set(confAddr, confSpaceBuf.getPos());
				var atomMolisAddr = confSpaceBuf.place(arrayStruct);
				arrayStruct.size.set(atomMolisAddr, conf.coords.size);
				arrayStruct.things_ptr.set(atomMolisAddr, 0);
				for (int i=0; i<conf.coords.size; i++) {
					confSpaceBuf.int32(conf.atomMolInfoIndices[i]);
				}
				confSpaceBuf.skipToAlignment(WidestGpuLoad);

				// write the motions
				confStruct.motions_offset.set(confAddr, confSpaceBuf.getPos());
				var confMotionOffsetsAddr = confSpaceBuf.place(confMotionOffsets, conf.motions.length, WidestGpuLoad);
				for (int i=0; i<conf.motions.length; i++) {
					var motion = conf.motions[i];
					confMotionOffsets.set(confMotionOffsetsAddr, i, confSpaceBuf.getPos());
					if (motion instanceof DihedralAngle.Description) {
						writeDihedral((DihedralAngle.Description)motion, pos.index, confSpaceBuf);
					} else {
						throw new UnsupportedOperationException(motion.getClass().getName());
					}
				}

				assert (confSpaceBuf.isAligned(WidestGpuLoad));
			}
		}

		assert (confSpaceBuf.getPos() == confSpaceSize + positionsSize);
		assert (confSpaceBuf.isAligned(WidestGpuLoad));

		// write the static atom coords
		confSpaceStruct.static_atom_coords_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		var atomCoordsAddr = confSpaceBuf.place(arrayStruct);
		arrayStruct.size.set(atomCoordsAddr, confSpace.staticCoords.size);
		arrayStruct.things_ptr.set(atomCoordsAddr, 0);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			var addr = confSpaceBuf.place(real3Struct);
			real3Struct.x.set(addr, confSpace.staticCoords.x(i));
			real3Struct.y.set(addr, confSpace.staticCoords.y(i));
			real3Struct.z.set(addr, confSpace.staticCoords.z(i));
		}
		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		// write the static atom molecule indices
		confSpaceStruct.static_atom_molis_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		var atomMolisAddr = confSpaceBuf.place(arrayStruct);
		arrayStruct.size.set(atomMolisAddr, confSpace.staticCoords.size);
		arrayStruct.things_ptr.set(atomMolisAddr, 0);
		for (int i=0; i<confSpace.staticCoords.size; i++) {
			confSpaceBuf.int32(confSpace.staticMolInfoIndices[i]);
		}
		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		assert (confSpaceBuf.getPos() == confSpaceSize + positionsSize + staticCoordsSize);
		assert (confSpaceBuf.isAligned(WidestGpuLoad));

		// write the forcefield params
		confSpaceStruct.params_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		forcefieldsImpl.writeParams(confSpaceBuf);

		// write the pos pairs
		confSpaceStruct.pos_pairs_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		var posPairOffsetsAddr = confSpaceBuf.place(posPairOffsets, numPosPairs, WidestGpuLoad);

		assert (confSpaceBuf.isAligned(WidestGpuLoad));

		// write the static-static pair
		posPairOffsets.set(posPairOffsetsAddr, 0, confSpaceBuf.getPos());
		forcefieldsImpl.writeStaticStatic(confSpaceBuf);

		// write the static-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(posPairOffsetsAddr, 1 + posi1, confSpaceBuf.getPos());
			var fragOffsetsAddr = confSpaceBuf.place(fragOffsets, confSpace.numFrag(posi1), WidestGpuLoad);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragOffsetsAddr, fragi1, confSpaceBuf.getPos());
				forcefieldsImpl.writeStaticPos(posi1, fragi1, confSpaceBuf);
			}
		}

		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		// write the pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			posPairOffsets.set(posPairOffsetsAddr, 1 + confSpace.positions.length + posi1, confSpaceBuf.getPos());
			var fragOffsetsAddr = confSpaceBuf.place(fragOffsets, confSpace.numFrag(posi1), WidestGpuLoad);
			for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
				fragOffsets.set(fragOffsetsAddr, fragi1, confSpaceBuf.getPos());
				forcefieldsImpl.writePos(posi1, fragi1, confSpaceBuf);
			}
		}

		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		// write the pos-pos pairs
		for (int posi1=0; posi1<confSpace.positions.length; posi1++) {
			for (int posi2=0; posi2<posi1; posi2++) {
				posPairOffsets.set(posPairOffsetsAddr, 1 + 2*confSpace.positions.length + posi1*(posi1 - 1)/2 + posi2, confSpaceBuf.getPos());
				var fragOffsetsAddr = confSpaceBuf.place(fragOffsets, confSpace.numFrag(posi1)*confSpace.numFrag(posi2), WidestGpuLoad);
				for (int fragi1=0; fragi1<confSpace.numFrag(posi1); fragi1++) {
					for (int fragi2=0; fragi2<confSpace.numFrag(posi2); fragi2++) {
						fragOffsets.set(fragOffsetsAddr, fragi1*confSpace.numFrag(posi2) + fragi2, confSpaceBuf.getPos());
						forcefieldsImpl.writePosPos(posi1, fragi1, posi2, fragi2, confSpaceBuf);
					}
				}
			}
		}

		assert (confSpaceBuf.getPos() == confSpaceSize + positionsSize + staticCoordsSize + forcefieldSize);
		confSpaceBuf.skipToAlignment(WidestGpuLoad);

		// write the molecule motions
		confSpaceStruct.molecule_motions_offset.set(confSpaceAddr, confSpaceBuf.getPos());
		var molMotionOffsetsAddr = confSpaceBuf.place(molMotionOffsets, numMolMotions);
		int molMotionIndex = 0;
		for (int moli=0; moli<confSpace.molInfos.length; moli++) {
			var molInfo = confSpace.molInfos[moli];
			for (var motion : molInfo.motions) {
				molMotionOffsets.set(molMotionOffsetsAddr, molMotionIndex++, confSpaceBuf.getPos());
				if (motion instanceof DihedralAngle.Description) {
					writeDihedral((DihedralAngle.Description)motion, PosInter.StaticPos, confSpaceBuf);
				} else if (motion instanceof TranslationRotation.Description) {
					writeTranslationRotation((TranslationRotation.Description)motion, moli, confSpaceBuf);
				} else {
					throw new UnsupportedOperationException(motion.getClass().getName());
				}
			}
		}

		assert (confSpaceBuf.getPos() == confSpaceSize + positionsSize + staticCoordsSize + forcefieldSize + molMotionsSize);

		// make sure we used the whole buffer
		assert (confSpaceBuf.getPos() == bufSize) : String.format("%d bytes leftover", bufSize - confSpaceBuf.getPos());

		// upload the conf space to the GPU for each device
		var confSpaceBuf = this.confSpaceBuf.asByteBuffer();
		for (var it : gpuStreams) {
			Pointer p = switch (precision) {
				case Float32 -> NativeLib.alloc_conf_space_f32(it.gpuInfo.id, confSpaceBuf);
				case Float64 -> NativeLib.alloc_conf_space_f64(it.gpuInfo.id, confSpaceBuf);
			};
			pConfSpace.put(it.gpuInfo, p);
		}

		// make the conf space sizes struct
		confSpaceSizesMem = new MemoryBuffer(confSpaceSizesStruct.bytes());
		confSpaceSizesStruct.num_pos.set(confSpaceSizesMem, confSpace.numPos());
		confSpaceSizesStruct.max_num_inters.set(confSpaceSizesMem, numPosPairs);
		confSpaceSizesStruct.num_atoms.set(confSpaceSizesMem, confSpace.maxNumConfAtoms);
		confSpaceSizesStruct.max_num_dofs.set(confSpaceSizesMem, confSpace.maxNumDofs);
		confSpaceSizesStruct.num_mol_motions.set(confSpaceSizesMem, numMolMotions);
		confSpaceSizesBuf = confSpaceSizesMem.asByteBuffer();

		// keep track of how many streams are left to give out
		class StreamsLeft {
			final GpuInfo gpuInfo;
			int numStreams;
			StreamsLeft(GpuStreams gpuStreams) {
				this.gpuInfo = gpuStreams.gpuInfo;
				this.numStreams = gpuStreams.numStreams;
			}
		}
		var gpuStreamsLeft = gpuStreams.stream()
			.map(StreamsLeft::new)
			.collect(Collectors.toList());

		// how big should the stream transfer buffers be?
		long hostBytes = switch (precision) {
			case Float32 -> NativeLib.minimize_batch_bufsize_host_f32(confSpaceSizesMem.asByteBuffer(), maxBatchSize);
			case Float64 -> NativeLib.minimize_batch_bufsize_host_f64(confSpaceSizesMem.asByteBuffer(), maxBatchSize);
		};
		long deviceBytes = switch (precision) {
			case Float32 -> NativeLib.minimize_batch_bufsize_device_f32(confSpaceSizesMem.asByteBuffer(), maxBatchSize);
			case Float64 -> NativeLib.minimize_batch_bufsize_device_f64(confSpaceSizesMem.asByteBuffer(), maxBatchSize);
		};

		// layout the gpu streams in a round-robin order
		int gpui = 0;
		while (true) {

			// find the next device with a stream left, if any
			StreamsLeft hasStreams = null;
			for (int i=0; i<gpuStreams.size() && hasStreams == null; i++) {
				var it = gpuStreamsLeft.get((gpui + i) % gpuStreams.size());
				if (it.numStreams > 0) {
					hasStreams = it;
				}
			}

			// if no streams left to distribute, we're done here
			if (hasStreams == null) {
				break;
			}

			// allocate a stream
			hasStreams.numStreams -= 1;
			streams.add(new Stream(hasStreams.gpuInfo.id, pConfSpace.get(hasStreams.gpuInfo), hostBytes, deviceBytes));

			// try the next gpu first next time
			gpui += 1;
		}

		// init the stream queue
		streamQueue = new ArrayBlockingQueue<>(streams.size());
		streamQueue.addAll(streams);
	}

	private void writeDihedral(DihedralAngle.Description desc, int posi, MemoryBuffer buf) {

		var dihedralAddr = buf.place(dihedralStruct);
		dihedralStruct.id.set(dihedralAddr, dihedralId);
		dihedralStruct.min_radians.set(dihedralAddr, Math.toRadians(desc.minDegrees));
		dihedralStruct.max_radians.set(dihedralAddr, Math.toRadians(desc.maxDegrees));
		dihedralStruct.a_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.a));
		dihedralStruct.b_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.b));
		dihedralStruct.c_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.c));
		dihedralStruct.d_index.set(dihedralAddr, desc.getAtomIndex(confSpace, posi, desc.d));
		dihedralStruct.num_rotated.set(dihedralAddr, desc.rotated.length);
		dihedralStruct.modified_posi.set(dihedralAddr, posi);

		var rotatedIndices = int32array();
		var indicesAddr = buf.place(rotatedIndices, desc.rotated.length, WidestGpuLoad);
		for (int i=0; i<desc.rotated.length; i++) {
			rotatedIndices.set(indicesAddr, i, desc.getAtomIndex(confSpace, posi, desc.rotated[i]));
		}

		buf.skipToAlignment(WidestGpuLoad);
	}

	private void writeTranslationRotation(TranslationRotation.Description desc, int moli, MemoryBuffer buf) {

		var transrotAddr = buf.place(transRotStruct);
		transRotStruct.id.set(transrotAddr, transRotId);
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

	@Override
	public ConfSpace confSpace() {
		return confSpace;
	}

	private static class Stream implements AutoCloseable {

		public final int device;
		public final Pointer pConfSpace;
		public final Pointer stream;

		public Stream(int device, Pointer pConfSpace, long hostBytes, long deviceBytes) {
			this.device = device;
			this.pConfSpace = pConfSpace;
			this.stream = NativeLib.alloc_stream(device, hostBytes, deviceBytes);
		}

		@Override
		public void close() {
			NativeLib.free_stream(device, stream);
		}
	}

	private final List<Stream> streams = new ArrayList<>();
	private final BlockingQueue<Stream> streamQueue;

	public class CheckedOutStream implements AutoCloseable {

		public final Stream stream;

		private CheckedOutStream(Stream stream) {
			this.stream = stream;
		}

		@Override
		public void close() {
			if (stream != null) {
				streamQueue.add(stream);
			}
		}
	}

	/**
	 * Check out a GPU stream.
	 * If none are immediately available, wait until one is.
	 */
	public CheckedOutStream checkoutStream() {

		Stream stream = null;
		try {
			while (stream == null) {
				stream = streamQueue.poll(100, TimeUnit.MILLISECONDS);
			}
		} catch (InterruptedException ex) {
			throw new RuntimeException(ex);
		}

		return new CheckedOutStream(stream);
	}

	/**
	 * Try to check out a GPU stream.
	 * If none are immediately available, the stream field will be null.
	 */
	public CheckedOutStream tryCheckoutStream() {

		// try to get a stream
		return new CheckedOutStream(streamQueue.poll());
	}

	public int numStreams() {
		return streams.size();
	}

	@Override
	public void close() {
		for (var it : gpuStreams) {
			NativeLib.free_conf_space(it.gpuInfo.id, pConfSpace.get(it.gpuInfo));
		}
		confSpaceSizesMem.close();
		for (Stream stream : streams) {
			stream.close();
		}
		confSpaceBuf.close();
	}

	public AssignedCoords assign(int[] conf) {
		try (var confMem = makeConf(conf)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				try (var co = checkoutStream()) {
					switch (precision) {
						case Float32 -> NativeLib.assign_f32(co.stream.device, co.stream.stream, co.stream.pConfSpace, confMem.asByteBuffer(), coordsMem.asByteBuffer());
						case Float64 -> NativeLib.assign_f64(co.stream.device, co.stream.stream, co.stream.pConfSpace, confMem.asByteBuffer(), coordsMem.asByteBuffer());
					}
				}
				return makeCoords(coordsMem, conf);
			}
		}
	}

	private MemoryBuffer makeConf(int[] conf) {

		MemoryBuffer mem = makeArray(confSpace.positions.length, Int32.bytes);
		var array = int32array();
		var addr = getArrayAddress(mem);
		for (int posi=0; posi<confSpace.positions.length; posi++) {
			array.set(addr, posi, conf[posi]);
		}

		return mem;
	}

	private AssignedCoords makeCoords(MemoryBuffer mem, int[] assignments) {

		AssignedCoords coords = new AssignedCoords(confSpace, assignments);

		// copy the coords from the native memory
		var arrayAddr = getArrayAddress(mem);
		for (int i=0; i<confSpace.maxNumConfAtoms; i++) {
			var addr = arrayAddr.sliceFrom(i*real3Struct.bytes());
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
	public EnergiedCoords calc(int[] conf, List<PosInter> inters) {
		try (var confMem = makeConf(conf)) {
			try (var intersMem = makeIntersMem(inters)) {
				try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
					try (var co = checkoutStream()) {
						var energy = forcefieldsImpl.calc(co.stream, confMem.asByteBuffer(), intersMem.asByteBuffer(), coordsMem.asByteBuffer(), confSpace.maxNumConfAtoms);
						return new EnergiedCoords(
							makeCoords(coordsMem, conf),
							energy
						);
					}
				}
			}
		}
	}

	@Override
	public double calcEnergy(int[] conf, List<PosInter> inters) {
		try (var confMem = makeConf(conf)) {
			try (var intersMem = makeIntersMem(inters)) {
				try (var co = checkoutStream()) {
					return forcefieldsImpl.calc(co.stream, confMem.asByteBuffer(), intersMem.asByteBuffer(), null, confSpace.maxNumConfAtoms);
				}
			}
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
		List<MinimizationJob> jobs = Collections.singletonList(new MinimizationJob(conf, inters));
		try (var jobsMem = makeMinimizationJobsMem(jobs)) {
			try (var coordsMem = makeArray(confSpace.maxNumConfAtoms, real3Struct.bytes())) {
				try (var dofsMem = makeArray(confSpace.maxNumDofs, precision.bytes)) {
					try (var co = checkoutStream()) {
						double energy = forcefieldsImpl.minimize(co.stream, jobsMem.asByteBuffer(), coordsMem.asByteBuffer(), dofsMem.asByteBuffer());
						return new EnergiedCoords(
							makeCoords(coordsMem, conf),
							energy,
							makeDofs(dofsMem)
						);
					}
				}
			}
		}
	}

	@Override
	public double minimizeEnergy(int[] conf, List<PosInter> inters) {
		MinimizationJob job = new MinimizationJob(conf, inters);
		minimizeEnergies(Collections.singletonList(job));
		return job.energy;
	}

	@Override
	public void minimizeEnergies(List<MinimizationJob> jobs) {
		minimizeEnergies(jobs, null);
	}

	public void minimizeEnergies(List<MinimizationJob> jobs, CheckedOutStream co) {

		// check the batch size
		if (jobs.size() <= 0) {
			// nothing to do
			// NOTE: don't try to send an empty job list to the GPU minimizer, it will complain
			return;
		} else if (jobs.size() > maxBatchSize) {
			throw new IllegalArgumentException(String.format("too many jobs for batch: %d, max allocated is %d",
				jobs.size(), maxBatchSize
			));
		}

		try (var jobsMem = makeMinimizationJobsMem(jobs)) {

			if (dumpJobsTo != null) {
				forcefieldsImpl.dumpJobs(dumpJobsTo, jobsMem.asByteBuffer());
				throw new Error("Jobs dumped to: " + dumpJobsTo.getAbsolutePath());
			}

			try (var energiesMem = makeArray(jobs.size(), precision.bytes)) {

				// if we have a stream already, use that
				if (co != null) {
					forcefieldsImpl.minimizeBatch(co.stream, jobsMem.asByteBuffer(), energiesMem.asByteBuffer());
				} else {
					// otherwise, check one out now
					try (var co2 = checkoutStream()) {
						forcefieldsImpl.minimizeBatch(co2.stream, jobsMem.asByteBuffer(), energiesMem.asByteBuffer());
					}
				}

				// read the energies back into the jobs
				Real.Array energies = realarray(precision);
				var addr = getArrayAddress(energiesMem);
				for (int i=0; i<jobs.size(); i++) {
					jobs.get(i).energy = energies.get(addr, i);
				}
			}
		}
	}

	@Override
	public int maxBatchSize() {
		return (int)maxBatchSize;
	}

	private MemoryBuffer makeMinimizationJobsMem(List<MinimizationJob> jobs) {

		long memBytes =
			getArrayBytes(jobs.size(), Int64.bytes) // job offsets
			+ jobs.stream()
				.mapToLong(job ->
					getArrayBytes(confSpace.positions.length, Int32.bytes) // conf
					+ getArrayBytes(job.inters.size(), posInterStruct.bytes()) // inters
				)
				.sum();

		MemoryBuffer buf = new MemoryBuffer(memBytes);

		var addr = buf.place(arrayStruct);
		arrayStruct.size.set(addr, jobs.size());

		assert (buf.isAligned(WidestGpuLoad));

		var offsetsArray = int64array();
		var offsetsAddr = buf.place(offsetsArray, jobs.size(), WidestGpuLoad);

		assert (buf.isAligned(WidestGpuLoad));

		var confArray = int32array();
		for (int j=0; j<jobs.size(); j++) {
			var job = jobs.get(j);

			offsetsArray.set(offsetsAddr, j, buf.getPos());

			assert (buf.isAligned(WidestGpuLoad));

			// write the conf
			addr = buf.place(arrayStruct);
			arrayStruct.size.set(addr, job.conf.length);
			addr = buf.place(confArray, job.conf.length, WidestGpuLoad);
			for (int i=0; i<job.conf.length; i++) {
				confArray.set(addr, i, job.conf[i]);
			}

			assert (buf.isAligned(WidestGpuLoad));

			// write the inters
			addr = buf.place(arrayStruct);
			arrayStruct.size.set(addr, job.inters.size());
			for (var inter : job.inters) {
				addr = buf.place(posInterStruct);
				posInterStruct.posi1.set(addr, inter.posi1);
				posInterStruct.posi2.set(addr, inter.posi2);
				posInterStruct.weight.set(addr, inter.weight);
				posInterStruct.offset.set(addr, inter.offset);
			}
			buf.skipToAlignment(WidestGpuLoad);
		}

		// make sure we used the whole buffer
		assert (buf.getPos() == memBytes);

		return buf;
	}

	private MemoryBuffer makeIntersMem(List<PosInter> inters) {
		MemoryBuffer buf = makeArray(inters.size(), posInterStruct.bytes());
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

	private static long getArrayBytes(long size, long itemBytes) {
		return padToGpuAlignment(Int64.bytes*2 + size*itemBytes);
	}

	private static MemoryBuffer makeArray(long size, long itemBytes) {
		MemoryBuffer buf = new MemoryBuffer(getArrayBytes(size, itemBytes));
		buf.int64(size);
		buf.skip(8); // padding
		return buf.sliceFrom(0); // rewind the buffer, position is now at zero.
	}

	private static long getArraySize(MemoryBuffer mem) {
		return mem.getLong(0);
	}

	private static MemoryBuffer getArrayAddress(MemoryBuffer mem) {
		return mem.sliceFrom(Int64.bytes*2);
	}
}
