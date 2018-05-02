package edu.duke.cs.osprey.pruning;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.gmec.PruningSettings;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

import java.math.BigInteger;
import java.util.ArrayList;

public class BenchmarkSimpleDEE extends TestBase {

	public static void main(String[] args) {

		initDefaultEnvironment();

		//String flexibleResNumbers = "3 4 6";
		//String flexibleResNumbers = "3 4 6 7 8 9 10";
		String flexibleResNumbers = "3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20";
		//String flexibleResNumbers = "2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30";

		// make a new-school conf space
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : flexibleResNumbers.split(" ")) {
			strand.flexibility.get("A" + resNum).setLibraryRotamers(Strand.WildType, "ARG");
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			// calc an energy matrix
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			// run precomputed matrices DEE
			double pruningInterval = 10.0;
			double Ew = 0.0;
			String name = "matrices";
			PruningSettings pruningSettings = new PruningSettings();
			pruningSettings.algOption = 3;
			pruningSettings.useTriples = true;
			Stopwatch matsStopwatch = new Stopwatch().start();
			PrecomputedMatrices mats = new PrecomputedMatrices(
				pruningInterval, Ew, name, emat, confSpace,
				null, null, new EPICSettings(), new LUTESettings(), pruningSettings
			);
			mats.getPruneMat().prunePairsFromSingles();
			matsStopwatch.stop();
			System.out.println("precomputed matrices DEE finished in " + matsStopwatch.getTime(2));

			/* report on the conf space size
			BigInteger allConfs = new RCs(confSpace).getNumConformations();
			BigInteger unprunedLower = mats.getPruneMat().calcUnprunedConfsLowerBound();
			BigInteger unprunedUpper = mats.getPruneMat().calcUnprunedConfsUpperBound();
			BigInteger prunedLower = allConfs.subtract(unprunedUpper);
			BigInteger prunedUpper = allConfs.subtract(unprunedLower);
			double percentPrunedLower = 100.0*prunedLower.doubleValue()/allConfs.doubleValue();
			double percentPrunedUpper = 100.0*prunedUpper.doubleValue()/allConfs.doubleValue();
			log("Conformations defined by conformation space:                           %s", formatBig(allConfs));
			log("Conformations pruned (by singles and pairs, bounds):                   [%s,%s]", formatBig(prunedLower), formatBig(prunedUpper));
			log("Conformations remaining after pruning (by singles and pairs, bounds):  [%s,%s]", formatBig(unprunedLower), formatBig(unprunedUpper));
			log("Percent conformations pruned (by singles and pairs, bounds):           [%.6f,%.6f]", percentPrunedLower, percentPrunedUpper);
			*/

			// run SimpleDEE
			SimpleDEE.Runner runner = new SimpleDEE.Runner()
				.setThreshold(pruningSettings.stericThresh)
				.setSinglesGoldsteinDiffThreshold(pruningInterval)
				.setPairsGoldsteinDiffThreshold(pruningInterval)
				.setTriplesGoldsteinDiffThreshold(pruningInterval)
				.setParallelism(Parallelism.makeCpu(8))
				.setShowProgress(true);
			Stopwatch newStopwatch = new Stopwatch().start();
			PruningMatrix pmat = runner.run(confSpace, emat);
			newStopwatch.stop();
			System.out.print("simple DEE finished in " + newStopwatch.getTime(2));
			System.out.println(String.format(" speedup over precomputed matrices: %.2fx", matsStopwatch.getTimeMs()/newStopwatch.getTimeMs()));

			// TEMP
			//assertThat(pmat, is(mats.getPruneMat()));
			//System.out.println("pruning matrices are identical");
		}
	}
}
