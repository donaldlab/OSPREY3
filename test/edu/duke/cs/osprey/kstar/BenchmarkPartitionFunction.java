/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
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

import java.util.function.Supplier;


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

			RCs rcs = new RCs(confSpace);
			Supplier<ConfAStarTree> astarFactory = () -> new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();

			benchmarkPfunc(new SimplePartitionFunction(confEcalc, astarFactory.get(), rcs.getNumConformations()));
			benchmarkPfunc(new GradientDescentPfunc(confEcalc, astarFactory.get(), astarFactory.get(), rcs.getNumConformations()));
		}
	}

	private static void benchmarkPfunc(PartitionFunction pfunc) {

		// compute pfunc
		pfunc.init(0.05);
		Stopwatch sw = new Stopwatch().start();
		pfunc.compute();
		System.out.println(String.format("%-30s %s", pfunc.getClass().getSimpleName(), sw.stop().getTime(2)));
	}
}
