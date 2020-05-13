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
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.Structs;

import java.util.Arrays;
import java.util.List;

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

	private static void benchmark(String name, Benchmark[] bm, Benchmark base, int[] threadSizes, int numWarmups, int numRuns, Runnable task) {
		benchmark(name, bm, base, threadSizes, numWarmups, numRuns, task, () -> {});
	}

	private static void benchmark(String name, Benchmark[] bm, Benchmark base, int[] threadSizes, int numWarmups, int numRuns, Runnable task, Runnable cleanup) {
		log("%s:", name);
		for (int i=0; i<threadSizes.length; i++) {
			bm[i] = new Benchmark(threadSizes[i], numWarmups, numRuns, task);
			cleanup.run();
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

	private static void benchmarkMinimize(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

		log("Minimize:");

		// make the wild-type conformations
		int[] classicConf = classic.makeConfComplexWt();
		int[] compiledConf = compiled.makeConfComplexWt();

		// make interactions for the classic design case (ie, without the static-static contribution)
		ResidueInteractions classicInters = EnergyPartition.makeFragment(classic.complex, null, false, new RCTuple(classicConf));
		List<PosInter> compiledInters = PosInterDist.dynamic(compiled.complex);

		// NOTE: workbox has 12 cores, and 15 SMs
		// NOTE: jerrys have 48 cores, and 4x80 SMs
		int[] threadSizes = { 1, 3, 6, 12, 24, 48 };
		int[] streamSizes = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 320 };
		int numWarmups = 2;
		int numRuns = 10;

		Benchmark[] bmClassic = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiled = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledf32 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledf64 = new Benchmark[threadSizes.length];
		Benchmark[] bmCompiledCudaf32 = new Benchmark[streamSizes.length];
		Benchmark[] bmCompiledCudaf64 = new Benchmark[streamSizes.length];

		// classic
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(true)
			.build()) {

			benchmark("classic", bmClassic, null, threadSizes, numWarmups, numRuns, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
		}

		{ // compiled
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			benchmark("compiled", bmCompiled, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}

		{ // compiled native f32
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			benchmark("compiled f32", bmCompiledf32, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}
		{ // compiled native f64
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			benchmark("compiled f64", bmCompiledf64, bmClassic[0], threadSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
		}

		{ // compiled CUDA f32
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			benchmark("compiled CUDA f32", bmCompiledCudaf32, bmClassic[0], streamSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			}, ecalc::recycleStreams);
		}

		{ // compiled CUDA f64
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			benchmark("compiled CUDA f64", bmCompiledCudaf64, bmClassic[0], streamSizes, numWarmups, numRuns, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			}, ecalc::recycleStreams);
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

		for (Structs.Precision precision : Structs.Precision.values()) {
		//{ Structs.Precision precision = Structs.Precision.Float32;
		//{ Structs.Precision precision = Structs.Precision.Float64;
			try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {

				// use all the interactions
				List<PosInter> inters = PosInterDist.all(confSpace);

				int[] conf = new int[] { 0, 0, 0, 0, 0, 0, 0 };

				log("precision = %s", precision);

				log("energy = %f", confEcalc.calcEnergy(conf, inters));
				log("   exp = %f", 2199.44093411);

				log("energy = %f", confEcalc.minimizeEnergy(conf, inters));
				log("   exp = %f", -1359.27208010);
			}
		}
	}
}
