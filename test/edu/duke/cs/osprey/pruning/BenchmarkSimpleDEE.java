package edu.duke.cs.osprey.pruning;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tupexp.LUTESettings;

import java.util.ArrayList;

public class BenchmarkSimpleDEE extends TestBase {

	public static void main(String[] args) {

		initDefaultEnvironment();

		//String flexibleResNumbers = "3 4 6";
		String flexibleResNumbers = "3 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20";
		//String flexibleResNumbers = "2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30";

		// make a new-school conf space
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		for (String resNum : flexibleResNumbers.split(" ")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType);
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.use((ecalc) -> {

				// calc an energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcEnergyMatrix();

				// get an old-school search problem so we can use the original DEE code
				TestBase.ResidueFlexibility resFlex = new TestBase.ResidueFlexibility();
				resFlex.addFlexible(flexibleResNumbers);
				resFlex.sortPositions(); // to match the SimpleConfSpace
				boolean addWt = false;
				boolean useEpic = false;
				boolean useTupleExpansion = false;
				boolean useEllipses = false;
				boolean useERef = false;
				boolean addResEntropy = false;
				boolean addWtRots = false;
				boolean doMinimize = false;
				ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
				ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
				SearchProblem search = new SearchProblem(
					"test", "examples/1CC8/1CC8.ss.pdb",
					resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
					new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
					false, new ArrayList<>()
				);

				search.emat = emat;

				// setup old-school DEE (steric singles+pairs and goldstein singles+pairs)
				double pruningInterval = 100.0; // I0 + Ew
				boolean typeDep = false; // default = false
				double boundsThreshold = 0; // default = 100, but never used in DEE as far as I can tell
				int algOption = 3; // default = 1
				boolean useFlags = true; // default = true
				boolean useTriples = false; // default = false
				boolean preDACS = false; // default = false
				double stericThreshold = 100; // default = 100
				search.pruneMat = new PruningMatrix(search.confSpace, 0);
				PruningControl pruning = new PruningControl(search,
					pruningInterval, typeDep, boundsThreshold, algOption, useFlags,
					useTriples, preDACS, useEpic, useTupleExpansion, stericThreshold
				);

				// benchmark it
				Stopwatch oldStopwatch = new Stopwatch().start();
				pruning.prune();
				oldStopwatch.stop();
				System.out.println("old-school DEE finished in " + oldStopwatch.getTime(2));

				// update pruned pairs with singles info, so the old pmat can match the new pmat
				search.pruneMat.prunePairsFromSingles();

				// setup new-school DEE
				SimpleDEE.Runner runner = new SimpleDEE.Runner()
					.setThreshold(100.0)
					.setGoldsteinDiffThreshold(100.0)
					.setShowProgress(true);

				// benchmark it
				Stopwatch newStopwatch = new Stopwatch().start();
				PruningMatrix pmat = runner.run(confSpace, emat);
				newStopwatch.stop();
				System.out.print("new-school DEE finished in " + newStopwatch.getTime(2));
				System.out.println(String.format(" speedup: %.2fx", oldStopwatch.getTimeMs()/newStopwatch.getTimeMs()));

				assertThat(pmat, is(search.pruneMat));
				System.out.println("pruning matrices are identical");
			});
	}
}
