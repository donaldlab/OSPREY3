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

import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkEnergies {

	public static void main(String[] args) {

		// make the conf spaces
		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		//benchmarkEcalc(classic, compiled);
		//benchmarkEmatCpu(classic, compiled);

		nativeLab(compiled);

		// TODO: with static-static energies on compiled ecalcs?
		// TODO: pfuncs
		// TODO: GPUs
	}

	private static void benchmarkEcalc(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

		// make the wild-type conformations
		int[] classicConf = classic.makeConfComplexWt();
		int[] compiledConf = compiled.makeConfComplexWt();

		// make interactions for the classic design case (ie, without the static-static contribution)
		ResidueInteractions classicInters = EnergyPartition.makeFragment(classic.complex, null, false, new RCTuple(classicConf));
		List<PosInter> compiledInters = PosInterDist.dynamic(compiled.complex);

		log("Rigid energy:");
		Benchmark bmClassicRigid;
		Benchmark bmCompiledRigid;
		Benchmark bmCompiledf32Rigid;
		Benchmark bmCompiledf64Rigid;
		Benchmark bmCompiledCudaf32Rigid;
		Benchmark bmCompiledCudaf64Rigid;

		int numWarmupsRigid = 100;
		int numRunsRigid = 4000;

		// benchmark classic rigid energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(false)
			.build()) {

			bmClassicRigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
			log("\t%20s: %s", "classic", bmClassicRigid.toString());
		}

		{ // benchmark compiled rigid energies
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			bmCompiledRigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled", bmCompiledRigid.toString(bmClassicRigid));
		}

		{ // benchmark compiled native f32 rigid energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledf32Rigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled f32", bmCompiledf32Rigid.toString(bmClassicRigid));
		}
		{ // benchmark compiled native f64 rigid energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledf64Rigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled f64", bmCompiledf64Rigid.toString(bmClassicRigid));
		}

		{ // benchmark compiled gpu f32 rigid energies
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledCudaf32Rigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled CUDA f32", bmCompiledCudaf32Rigid.toString(bmClassicRigid));
		}
		{ // benchmark compiled gpu f64 rigid energies
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledCudaf64Rigid = new Benchmark(numWarmupsRigid, numRunsRigid, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled CUDA f64", bmCompiledCudaf64Rigid.toString(bmClassicRigid));
		}


		log("Minimized energy:");
		Benchmark bmClassicMinimized;
		Benchmark bmCompiledMinimized;
		Benchmark bmCompiledf32Minimized;
		Benchmark bmCompiledf64Minimized;
		Benchmark bmCompiledCudaf32Minimized;
		Benchmark bmCompiledCudaf64Minimized;

		int numWarmupsMinimized = 2;
		int numRunsMinimized = 40;

		// benchmark classic minimized energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(true)
			.build()) {

			bmClassicMinimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
			log("\t%20s: %s", "classic", bmClassicMinimized.toString());
		}

		{ // benchmark compiled minimized energies
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			bmCompiledMinimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled", bmCompiledMinimized.toString(bmClassicMinimized));
		}

		{ // benchmark compiled native f32 minimized energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledf32Minimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled f32", bmCompiledf32Minimized.toString(bmClassicMinimized));
		}
		{ // benchmark compiled native f64 minimized energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledf64Minimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled f64", bmCompiledf64Minimized.toString(bmClassicMinimized));
		}

		{ // benchmark compiled CUDA f32 minimized energies
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledCudaf32Minimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled CUDA f32", bmCompiledCudaf32Minimized.toString(bmClassicMinimized));
		}
		{ // benchmark compiled CUDA f64 minimized energies
			CudaConfEnergyCalculator ecalc = new CudaConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledCudaf64Minimized = new Benchmark(numWarmupsMinimized, numRunsMinimized, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%20s: %s", "compiled CUDA f64", bmCompiledCudaf64Minimized.toString(bmClassicMinimized));
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

		ConfSpace confSpace = compiled.complex;

		for (Structs.Precision precision : Structs.Precision.values()) {
		//{ Structs.Precision precision = Structs.Precision.Float32;
		//{ Structs.Precision precision = Structs.Precision.Float64;
			try (var confEcalc = new CudaConfEnergyCalculator(confSpace, precision)) {

				// use all the interactions
				List<PosInter> inters = PosInterDist.all(confSpace);

				int[] conf = new int[] { 0, 0, 0, 0, 0, 0, 0 };

				log("precision = %s", precision);

				// TEMP
				log("energy = %f", confEcalc.calcEnergy(conf, inters));
				log("   exp = %f", 2199.44093411);

				log("energy = %f", confEcalc.minimizeEnergy(conf, inters));
				log("   exp = %f", -1359.27208010);
			}
		}
	}
}
