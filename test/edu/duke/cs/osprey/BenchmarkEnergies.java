package edu.duke.cs.osprey;


import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.confspace.compiled.TestConfSpace;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.TaskExecutor;

import java.util.List;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkEnergies {

	public static void main(String[] args) {

		// make the conf spaces
		TestConfSpace.AffinityClassic classic = TestConfSpace.Design2RL0Interface7Mut.makeClassic();
		TestConfSpace.AffinityCompiled compiled = TestConfSpace.Design2RL0Interface7Mut.makeCompiled();

		benchmarkEcalcCpu(classic, compiled);

		// TODO: energy matrices
		// TODO: pfuncs
		// TODO: GPUs
	}

	private static void benchmarkEcalcCpu(TestConfSpace.AffinityClassic classic, TestConfSpace.AffinityCompiled compiled) {

		// make the wild-type conformations
		int[] classicConf = classic.makeConfComplexWt();
		int[] compiledConf = compiled.makeConfComplexWt();

		// make interactions for the classic design case (ie, without the static-static contribution)
		ResidueInteractions classicInters = EnergyPartition.makeFragment(classic.complex, null, false, new RCTuple(classicConf));
		List<PosInter> compiledInters = PosInterDist.dynamic(compiled.complex, null, compiledConf);

		log("Rigid energy:");
		Benchmark bmClassicRigid;
		Benchmark bmCompiledRigid;

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

		// benchmark compiled rigid energies
		try (TaskExecutor tasks = new TaskExecutor()) {
			ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex, tasks);

			bmCompiledRigid = new Benchmark(100, 5000, () -> {
				ecalc.calcEnergy(compiledConf, compiledInters);
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

			bmClassicMinimized = new Benchmark(5, 80, () -> {
				ParametricMolecule pmol = classic.complex.makeMolecule(classicConf);
				ecalc.calcEnergy(pmol, classicInters);
			});
			log("\t%10s: %s", "classic", bmClassicMinimized.toString());
		}

		// benchmark compiled minimized energies
		try (TaskExecutor tasks = new TaskExecutor()) {
			ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(compiled.complex, tasks);

			bmCompiledMinimized = new Benchmark(5, 80, () -> {
				ecalc.minimizeEnergy(compiledConf, compiledInters);
			});
			log("\t%10s: %s", "compiled", bmCompiledMinimized.toString(bmClassicMinimized));
		}
	}
}
