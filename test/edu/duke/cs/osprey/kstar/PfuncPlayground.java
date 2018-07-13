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
			try (ConfDB confdb = new ConfDB(confEcalc.confSpace, confdbFile)) {
				ConfDB.ConfTable table = confdb.new ConfTable("pfunctest");

				estimate(confEcalc, emat, surf, table);
			}

			surf.writeTraces(new File("pfunctest.trace.vtk"));
		}
	}

	public static void estimate(ConfEnergyCalculator confEcalc, EnergyMatrix emat, PfuncSurface surf, ConfDB.ConfTable table) {

		final double epsilon = 0.01;

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confEcalc.confSpace)
			.setTraditional()
			.build();

		GradientDescentPfunc pfunc = new GradientDescentPfunc(confEcalc);
		pfunc.init(astar, astar.getNumConformations(), epsilon);
		pfunc.traceTo(surf);
		pfunc.compute(surf.numEnergies);
	}
}
