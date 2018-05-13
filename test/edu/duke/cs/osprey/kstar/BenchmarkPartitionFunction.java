/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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

			benchmarkPfunc(confSpace, emat, new SimplePartitionFunction(confEcalc));
			benchmarkPfunc(confSpace, emat, new GradientDescentPfunc(confEcalc));
		}
	}

	private static void benchmarkPfunc(SimpleConfSpace confSpace, EnergyMatrix emat, PartitionFunction pfunc) {

		// make the A* search
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build();

		// compute pfunc
		final double targetEpsilon = 0.05;
		pfunc.init(astar, astar.getNumConformations(), targetEpsilon);
		Stopwatch sw = new Stopwatch().start();
		pfunc.compute();
		System.out.println(String.format("%-30s %s", pfunc.getClass().getSimpleName(), sw.stop().getTime(2)));
	}
}
