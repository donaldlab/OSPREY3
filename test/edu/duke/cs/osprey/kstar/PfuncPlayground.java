package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PfuncSurface;
import edu.duke.cs.osprey.parallelism.Parallelism;

import java.io.File;


public class PfuncPlayground {

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

			final int scoreBatch = 200;
			final int numEnergies = 1000;
			final int numScoreBatches = 1000;

			PfuncSurface surf = new PfuncSurface(scoreBatch, numScoreBatches, numEnergies);
			//surf.sample(confEcalc, emat);
			//surf.write(new File("pfunctest.vtk"));

			final File confdbFile = new File("pfunctext.conf.db");
			new ConfDB(confEcalc.confSpace, confdbFile).use((confdb) -> {
				ConfDB.ConfTable table = confdb.new ConfTable("pfunctest");

				estimate(confEcalc, emat, surf, table);
			});

			surf.writeTraces(new File("pfunctest.trace.vtk"));
		}
	}

	public static void estimate(ConfEnergyCalculator confEcalc, EnergyMatrix emat, PfuncSurface surf, ConfDB.ConfTable table) {

		final double epsilon = 0.01;

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confEcalc.confSpace)
			.setTraditional()
			.build();

		GradientDescentPfunc pfunc = new GradientDescentPfunc(astar, confEcalc);
		pfunc.init(epsilon);
		pfunc.traceTo(surf);
		pfunc.compute(surf.numEnergies);
	}
}
