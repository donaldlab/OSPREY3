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

package edu.duke.cs.osprey.energy.forcefield;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrixCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.Arrays;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class BenchmarkApproximatedForcefield {

	public static void main(String[] args) {

		// make a conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A16", "A17", "A18", "A19", "A20", "A21")) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		Function<ConfEnergyCalculator,EnergyMatrix> calcEmat = (confEcalc) ->
			new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setTripleCorrectionThreshold(10.0)
				.build()
				.calcEnergyMatrix();

		// get an energy calculator
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(1))
			.build()
		) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				//.setEnergyPartition(EnergyPartition.Traditional)
				.setEnergyPartition(EnergyPartition.AllOnPairs)
				.build();

			// calc the approximator matrix
			ApproximatorMatrix amat = new ApproximatorMatrixCalculator(confEcalc)
				.setNumSamplesPerDoF(9)
				.calc();

			ConfEnergyCalculator confEcalcApprox = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(confEcalc.epart)
				.setApproximatorMatrix(amat)
				.setApproximationErrorBudget(1e-2)
				.build();

			// calc energy matrices

			log("\ncalculating real emat...");
			Stopwatch realStopwatch = new Stopwatch().start();
			EnergyMatrix emat = calcEmat.apply(confEcalc);
			log("\tdone in %s", realStopwatch.stop().getTime(2));

			log("\ncalculating approx emat...");
			Stopwatch approxStopwatch = new Stopwatch().start();
			EnergyMatrix ematApprox = calcEmat.apply(confEcalcApprox);
			log("\tdone in %s", approxStopwatch.stop().getTime(2));
		}
	}
}
