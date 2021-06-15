package edu.duke.cs.osprey;


import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.compiled.*;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator.MinimizationJob;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.Structs;
import org.joml.Vector3d;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


public class BenchmarkEnergies {

	public static void main(String[] args) {

		// make the conf spaces
		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		benchmarkMinimize(classic, compiled);
		//benchmarkEmatCpu(classic, compiled);

		//nativeLab(compiled);

		// TODO: with static-static energies on compiled ecalcs?
		// TODO: pfuncs
	}

	private static void benchmarkThreads(String name, Benchmark[] bm, Benchmark base, int[] threadSizes, int numWarmups, int numRuns, Runnable task) {
		log("%s:", name);
		for (int i=0; i<threadSizes.length; i++) {

			// run the benchmark
			bm[i] = new Benchmark(threadSizes[i], numWarmups, numRuns, 1, task);

			// show the results
			logf("\t%2d threads: %s", threadSizes[i], bm[i].toString());
			if (base != null) {
				logf(", speedup over base %.2fx", bm[i].opsPerSecond/base.opsPerSecond);
			}
			if (i > 0) {
				logf(", speedup over single %.2fx", bm[i].opsPerSecond/bm[0].opsPerSecond);
			}
			logf("\n");
		}
	}

	private static <T extends AutoCloseable> void benchmarkGpus(String name, Benchmark[] bm, Benchmark base, int[] gpuSizes, int threadsPerGpu, int batchSize, int numWarmups, int numRuns, Function<Integer,T> init, Consumer<T> task) {
		log("%s:", name);
		for (int i=0; i<gpuSizes.length; i++) {

			int numGpus = gpuSizes[i];
			int numThreads = numGpus*threadsPerGpu;

			// do the init
			try (T ctx = init.apply(numGpus)) {

				// run the benchmark
				bm[i] = new Benchmark(numThreads, numWarmups, numRuns, batchSize, () -> task.accept(ctx));

			} catch (Exception ex) {
				throw new RuntimeException(ex);
			}

			// show the results
			logf("\t%1d gpus (%2d threads): %s", gpuSizes[i], numThreads, bm[i].toString());
			if (base != null) {
				logf(", speedup over base %.2fx", bm[i].opsPerSecond/base.opsPerSecond);
			}
			if (i > 0) {
				logf(", speedup over single %.2fx", bm[i].opsPerSecond/bm[0].opsPerSecond);
			}
			logf("\n");
		}
	}


	private static void benchmarkMinimize(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

		log("Minimize:");

		// make the wild-type conformations
		int[] classicConf = classic.makeConfWt(classic.complex);
		int[] compiledConf = compiled.makeConfWt(compiled.complex);

		// make interactions for the classic design case (ie, without the static-static contribution)
		ResidueInteractions classicInters = EnergyPartition.makeFragment(classic.complex, null, false, new RCTuple(classicConf));
		List<PosInter> compiledInters = PosInterDist.dynamic(compiled.complex, compiledConf);

		// NOTE: the jerrys have 48 cores, 4 GPUs
		int[] threadSizes = { 1, 3, 6, 12, 24, 48 };
		int[] gpuSizes = { 1, 2, 4 };
		int numWarmups = 2;
		int numRuns = 10;

		// assume all the GPUs are the same to get the best batch size
		var gpus = CudaConfEnergyCalculator.getGpusInfos();
		int batchSize = gpus.get(0).bestBatchSize();
		int threadsPerGpu = gpus.get(0).bestNumStreams();

		// set up the minimization jobs
		var jobs = new ArrayList<MinimizationJob>();
		for (int i=0; i<batchSize; i++) {
			jobs.add(new MinimizationJob(compiledConf, compiledInters));
		}

		Benchmark[] bmClassic = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiled = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledf32 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledf64 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledIntelf32 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledIntelf64 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledCudaf32 = new Benchmark[gpuSizes.length];
		Benchmark[] bmCompiledCudaf64 = new Benchmark[gpuSizes.length];

		// classic
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(true)
			.build()) {

			benchmarkThreads("classic", bmClassic, null, threadSizes, numWarmups, numRuns, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
		}

		{ // compiled
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			benchmarkThreads("compiled", bmCompiled, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}

		{ // compiled native f32
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			benchmarkThreads("compiled reference f32", bmCompiledf32, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}
		{ // compiled native f64
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			benchmarkThreads("compiled reference f64", bmCompiledf64, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}

		{ // compiled CUDA f32
			// NOTE: this benchmark is slower than it should be because the CCD minimizer gets unlucky on these inputs
			// if you tweak the step size for the line search a bit, f32 goes MUCH faster!
			benchmarkGpus("compiled CUDA f32", bmCompiledCudaf32, bmClassic[0], gpuSizes, threadsPerGpu, batchSize, numWarmups, numRuns,
				numGpus -> new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float32, gpus.subList(0, numGpus)),
				ecalc -> ecalc.minimizeEnergies(jobs)
			);
		}

		{ // compiled CUDA f64
			benchmarkGpus("compiled CUDA f64", bmCompiledCudaf64, bmClassic[0], gpuSizes, threadsPerGpu, batchSize, numWarmups, numRuns,
				numGpus -> new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float64, gpus.subList(0, numGpus)),
				ecalc -> ecalc.minimizeEnergies(jobs)
			);
		}
	}

	private static void benchmarkEmatCpu(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

		log("Rigid energy:");
		Benchmark bmClassicRigid;
		Benchmark bmCompiledRigid;

		// benchmark classic rigid energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(false)
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(classic.complex, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional)
				.build();
			SimplerEnergyMatrixCalculator ematCalc = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build();

			bmClassicRigid = new Benchmark(1, 4, () -> {
				ematCalc.calcEnergyMatrix();
			});
			log("\t%20s: %s", "classic", bmClassicRigid.toString());
		}

		// benchmark compiled rigid energies
		{
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(compiled.complex);
			EmatCalculator ematCalc = new EmatCalculator.Builder(confEcalc)
				.setIncludeStaticStatic(false)
				.setMinimize(false)
				.setPosInterDist(PosInterDist.DesmetEtAl1992)
				.build();

			bmCompiledRigid = new Benchmark(1, 4, () -> {
				ematCalc.calc();
			});
			log("\t%20s: %s", "compiled", bmCompiledRigid.toString(bmClassicRigid));
		}


		log("Minimized energy:");
		Benchmark bmClassicMinimized;
		Benchmark bmCompiledMinimized;

		// benchmark classic minimized energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(true)
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(classic.complex, ecalc)
				.setEnergyPartition(EnergyPartition.Traditional)
				.build();
			SimplerEnergyMatrixCalculator ematCalc = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build();

			bmClassicMinimized = new Benchmark(1, 4, () -> {
				ematCalc.calcEnergyMatrix();
			});
			log("\t%20s: %s", "classic", bmClassicMinimized.toString());
		}

		// benchmark compiled minimized energies
		{
			CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(compiled.complex);
			EmatCalculator ematCalc = new EmatCalculator.Builder(confEcalc)
				.setIncludeStaticStatic(false)
				.setMinimize(true)
				.setPosInterDist(PosInterDist.DesmetEtAl1992)
				.build();

			bmCompiledMinimized = new Benchmark(1, 4, () -> {
				ematCalc.calc();
			});
			log("\t%20s: %s", "compiled", bmCompiledMinimized.toString(bmClassicMinimized));
		}
	}

	private static void nativeLab(TestConfSpace.AffinityCompiled compiled) {

		List<CudaConfEnergyCalculator.GpuInfo> gpus = CudaConfEnergyCalculator.getGpusInfos();
		log("found %d GPUs", gpus.size());
		for (var gpu : gpus) {
			log(gpu.toString());
		}

		ConfSpace confSpace = compiled.complex;

		/* TEMP
		var gpuStreams = Collections.singletonList(
			new CudaConfEnergyCalculator.GpuStreams(CudaConfEnergyCalculator.getGpusInfos().get(0), 1)
		);
		int batchSize = 10;
		*/

		//for (Structs.Precision precision : Structs.Precision.values()) {
		//{ Structs.Precision precision = Structs.Precision.Float32;
		{ Structs.Precision precision = Structs.Precision.Float64;
			//try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision, gpuStreams, batchSize)) {
			try (var confEcalc = new NativeConfEnergyCalculator(confSpace, precision)) {

				int[] conf = new int[] { 0, 0, 0, 0, 0, 0, 0 };

				// use all the interactions
				List<PosInter> inters = PosInterDist.all(confSpace, conf);

				log("precision = %s", precision);

				//log("energy = %f", confEcalc.calcEnergy(conf, inters));
				//log("   exp = %f", 2199.44093411);

				//log("energy = %f", confEcalc.minimizeEnergy(conf, inters));

				/* TEMP
				var jobs = new ArrayList<MinimizationJob>(batchSize);
				for (int i=0; i<batchSize; i++) {
					jobs.add(new MinimizationJob(conf, inters));
				}
				confEcalc.minimizeEnergies(jobs);
				for (int i=0; i<batchSize; i++) {
					log("energy[%2d] = %f", i, jobs.get(i).energy);
				}
				*/

				var econf = confEcalc.minimize(conf, inters);
				log("energy = %f", econf.energy);

				Vector3d pos = new Vector3d();
				for (int i=0; i<10; i++) {
					econf.coords.coords.get(i, pos);
					log("\tpos[%2d] = %8.3f,%8.3f,%8.3f", i, pos.x, pos.y, pos.z);
				}
				for (int i=0; i<econf.dofValues.size(); i++) {
					log("\tdofs[%2d] = %12.6f", i, econf.dofValues.get(i));
				}

				//log("   exp = %f", -1359.27208010);
			}
		}
	}
}
