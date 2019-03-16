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
