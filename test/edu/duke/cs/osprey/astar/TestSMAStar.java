package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;


public class TestSMAStar {

	@Test
	public void tiny1CCC8() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		// use the minimum number of nodes, to make the harshest test for SMA*
		test(4, new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build()
		);
	}

	@Test
	public void small1CCC8() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		test(7, new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build()
		);
	}

	@Test
	public void medium1CCC8() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12")) {
			strand.flexibility.get(resNum).setLibraryRotamers("VAL");
		}

		// min number of nodes takes too long, so use a bit more
		test(1000, new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build()
		);
	}

	private static void test(int maxNumNodes, SimpleConfSpace confSpace) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {

			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();

			RCs rcs = new RCs(confSpace);

			// enumerate all the confs using A*
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.build();
			Stopwatch astarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> astarConfs = astar.nextConfs(Double.POSITIVE_INFINITY);
			astarStopwatch.stop();

			// enumerate all the confs using SMA*
			ConfAStarTree smastar = new ConfAStarTree.Builder(emat, rcs)
				.setTraditional()
				.setMaxNumNodes(maxNumNodes)
				.build();
			Stopwatch smastarStopwatch = new Stopwatch().start();
			List<ConfSearch.ScoredConf> smastarConfs = smastar.nextConfs(Double.POSITIVE_INFINITY);
			smastarStopwatch.stop();

			checkConfs(astarConfs, smastarConfs);
		}
	}

	private static void checkConfs(List<ConfSearch.ScoredConf> expectedConfs, List<ConfSearch.ScoredConf> observedConfs) {

		assertThat(observedConfs.size(), is(expectedConfs.size()));

		for (int i=0; i<expectedConfs.size(); i++) {
			ConfSearch.ScoredConf exp = expectedConfs.get(i);
			ConfSearch.ScoredConf obs = observedConfs.get(i);
			assertThat(obs.getAssignments(), is(exp.getAssignments()));
			assertThat(obs.getScore(), isAbsolutely(exp.getScore(), 1e-10));
		}
	}
}
