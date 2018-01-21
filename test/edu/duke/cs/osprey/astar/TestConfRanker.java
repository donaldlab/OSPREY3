package edu.duke.cs.osprey.astar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfRanker;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;


public class TestConfRanker {

	@Test
	public void tinyDiscrete1CC8() {

		// 8 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A7").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A9").setLibraryRotamers("GLY", "ALA");

		checkEveryConf(strand);
	}

	@Test
	public void small1CC8() {

		// 2268 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();

		checkEveryConf(strand);
	}

	private void checkEveryConf(Strand strand) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		// compute an emat
		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		// make the ranker
		ConfRanker ranker = new ConfRanker(confSpace, emat, (b) -> b.setTraditional());

		// check each conf has the correct rank, exhaustively
		int expectedNumConfs = 0;
		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build();
		while (true) {

			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}

			expectedNumConfs++;

			// roundoff error makes using the exact conf energies slightly unstable, so allow off-by-one errors
			int observedNumConfs = ranker.getNumConfsAtMost(conf.getScore()).intValueExact();
			assertThat(expectedNumConfs - observedNumConfs, lessThanOrEqualTo(1));

			// using a slight epsilon on both sides should be perfectly accurate though
			final double epsilon = 0.00001;
			assertThat(ranker.getNumConfsAtMost(conf.getScore() - epsilon).intValueExact(), is(expectedNumConfs - 1));
			assertThat(ranker.getNumConfsAtMost(conf.getScore() + epsilon).intValueExact(), is(expectedNumConfs));
		}
	}
}
