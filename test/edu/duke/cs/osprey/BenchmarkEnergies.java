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

		//benchmarkEcalcCpu(classic, compiled);
		//benchmarkEmatCpu(classic, compiled);

		nativeLab(compiled);

		// TODO: with static-static energies on compiled ecalcs?
		// TODO: pfuncs
		// TODO: GPUs
	}

	private static void benchmarkEcalcCpu(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

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

		// benchmark classic rigid energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(false)
			.build()) {

			bmClassicRigid = new Benchmark(100, 5000, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
			log("\t%10s: %s", "classic", bmClassicRigid.toString());
		}

		{ // benchmark compiled rigid energies
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			bmCompiledRigid = new Benchmark(100, 5000, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled", bmCompiledRigid.toString(bmClassicRigid));
		}

		{ // benchmark compiled native f32 rigid energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledf32Rigid = new Benchmark(100, 5000, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled f32", bmCompiledf32Rigid.toString(bmClassicRigid));
		}
		{ // benchmark compiled native f64 rigid energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledf64Rigid = new Benchmark(100, 5000, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled f64", bmCompiledf64Rigid.toString(bmClassicRigid));
		}


		log("Minimized energy:");
		Benchmark bmClassicMinimized;
		Benchmark bmCompiledMinimized;
		Benchmark bmCompiledf32Minimized;
		Benchmark bmCompiledf64Minimized;

		// benchmark classic minimized energies
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(classic.complex, new ForcefieldParams())
			.setIsMinimizing(true)
			.build()) {

			bmClassicMinimized = new Benchmark(5, 80, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
			log("\t%10s: %s", "classic", bmClassicMinimized.toString());
		}

		{ // benchmark compiled minimized energies
			CPUConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex);
			bmCompiledMinimized = new Benchmark(5, 80, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled", bmCompiledMinimized.toString(bmClassicMinimized));
		}

		{ // benchmark compiled native f32 minimized energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float32);
			bmCompiledf32Minimized = new Benchmark(5, 80, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled f32", bmCompiledf32Minimized.toString(bmClassicMinimized));
		}
		{ // benchmark compiled native f64 minimized energies
			NativeConfEnergyCalculator ecalc = new NativeConfEnergyCalculator(compiled.complex, Structs.Precision.Float64);
			bmCompiledf64Minimized = new Benchmark(5, 80, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled f64", bmCompiledf64Minimized.toString(bmClassicMinimized));
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
			log("\t%10s: %s", "classic", bmClassicRigid.toString());
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
			log("\t%10s: %s", "compiled", bmCompiledRigid.toString(bmClassicRigid));
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
			log("\t%10s: %s", "classic", bmClassicMinimized.toString());
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
			log("\t%10s: %s", "compiled", bmCompiledMinimized.toString(bmClassicMinimized));
		}
	}

	private static void nativeLab(TestConfSpace.AffinityCompiled compiled) {

		ConfSpace confSpace = compiled.complex;

		try (var confEcalc = new CudaConfEnergyCalculator(confSpace, Structs.Precision.Float64)) {

			// use all the interactions
			List<PosInter> inters = PosInterDist.all(confSpace);

			int[] conf = new int[] { 0, 0, 0, 0, 0, 0, 0 };

			double energy = confEcalc.calcEnergy(conf, inters);
			log("energy = %f", energy);
		}
	}
}
