package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;


public class BenchmarkPartitionFunction {

	public static void main(String[] args) {

		TestSimplePartitionFunction.TestInfo info = TestSimplePartitionFunction.make2RL0TestInfo();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(info.protein)
			.addStrand(info.ligand)
			.build();

		final Parallelism parallelism = Parallelism.make(8, 0, 0);
		final ForcefieldParams ffparams = new ForcefieldParams();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, ffparams)
			.setParallelism(parallelism)
			.build()
		) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc).build();
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.build()
				.calcEnergyMatrix();

			benchmarkPfunc(confSpace, emat, (astar) -> new SimplePartitionFunction(astar, confEcalc));
			benchmarkPfunc(confSpace, emat, (astar) -> new GradientDescentPfunc(astar, confEcalc));
		}
	}

	private static interface PfuncFactory {
		PartitionFunction make(ConfSearch astar);
	}

	private static void benchmarkPfunc(SimpleConfSpace confSpace, EnergyMatrix emat, PfuncFactory pfuncs) {

		// make the A* search
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build();

		// make the partition function
		PartitionFunction pfunc = pfuncs.make(astar);

		// compute pfunc
		final double targetEpsilon = 0.05;
		pfunc.init(targetEpsilon);
		Stopwatch sw = new Stopwatch().start();
		pfunc.compute();
		System.out.println(String.format("%-30s %s", pfunc.getClass().getSimpleName(), sw.stop().getTime(2)));
	}
}
