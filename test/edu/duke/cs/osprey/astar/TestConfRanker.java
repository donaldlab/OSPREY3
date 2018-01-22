package edu.duke.cs.osprey.astar;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ranking.*;
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

import java.math.BigInteger;
import java.util.Arrays;


public class TestConfRanker {

	public Strand makeTinyDiscrete1CC8() {

		// 8 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A7").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A9").setLibraryRotamers("GLY", "ALA");
		return strand;
	}
	@Test public void tinyDiscrete1CC8Sequential() { checkEveryConf(makeTinyDiscrete1CC8(), new SequentialOrderer()); }
	@Test public void tinyDiscrete1CC8DynamicOptimalPruning() { checkEveryConf(makeTinyDiscrete1CC8(), new DynamicOptimalPruningOrderer()); }
	@Test public void tinyDiscrete1CC8StaticOptimalPruning() { checkEveryConf(makeTinyDiscrete1CC8(), new StaticOptimalPruningOrderer()); }
	@Test public void tinyDiscrete1CC8DynamicHeuristicPruning() { checkEveryConf(makeTinyDiscrete1CC8(), new DynamicHeuristicPruningOrderer()); }

	public Strand makeSmall1CC8() {

		// 2268 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		return strand;
	}
	@Test public void small1CC8Sequential() { checkEveryConf(makeSmall1CC8(), new SequentialOrderer()); }
	@Test public void small1CC8DynamicOptimalPruning() { checkEveryConf(makeSmall1CC8(), new DynamicOptimalPruningOrderer()); }
	@Test public void small1CC8StaticOptimalPruning() { checkEveryConf(makeSmall1CC8(), new StaticOptimalPruningOrderer()); }
	@Test public void small1CC8DynamicHeuristicPruning() { checkEveryConf(makeSmall1CC8(), new DynamicHeuristicPruningOrderer()); }

	public Strand makeMedium1CC8() {

		// 3.62e6 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	@Test public void medium1CC8ZeroRankSequential() { assertThat(getZeroRank(makeMedium1CC8(), new SequentialOrderer()), is(new BigInteger("40306"))); }
	@Test public void medium1CC8ZeroRankDynamicOptimalPruning() { assertThat(getZeroRank(makeMedium1CC8(), new DynamicOptimalPruningOrderer()), is(new BigInteger("40306"))); }
	@Test public void medium1CC8ZeroRankStaticOptimalPruning() { assertThat(getZeroRank(makeMedium1CC8(), new StaticOptimalPruningOrderer()), is(new BigInteger("40306"))); }
	@Test public void medium1CC8ZeroRankDynamicHeuristicPruning() { assertThat(getZeroRank(makeMedium1CC8(), new DynamicHeuristicPruningOrderer()), is(new BigInteger("40306"))); }

	public Strand makeLarge1CC8() {

		// 3.62e6 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	@Test public void large1CC8ZeroRankSequential() { assertThat(getZeroRank(makeLarge1CC8(), new SequentialOrderer()), is(new BigInteger("1034629"))); }
	@Test public void large1CC8ZeroRankDynamicOptimalPruning() { assertThat(getZeroRank(makeLarge1CC8(), new DynamicOptimalPruningOrderer()), is(new BigInteger("1034629"))); }
	@Test public void large1CC8ZeroRankStaticOptimalPruning() { assertThat(getZeroRank(makeLarge1CC8(), new StaticOptimalPruningOrderer()), is(new BigInteger("1034629"))); }
	@Test public void large1CC8ZeroRankDynamicHeuristicPruning() { assertThat(getZeroRank(makeLarge1CC8(), new DynamicHeuristicPruningOrderer()), is(new BigInteger("1034629"))); }

	// NOTE: this one is big enough that's it's only useful for benchmarking
	public Strand makeHuge1CC8() {

		// 3.86e9 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	//@Test public void huge1CC8ZeroRankSequential() { assertThat(getZeroRank(makeHuge1CC8(), new SequentialOrderer()), is(new BigInteger("21039231"))); }
	//@Test public void huge1CC8ZeroRankDynamicOptimalPruning() { assertThat(getZeroRank(makeHuge1CC8(), new DynamicOptimalPruningOrderer()), is(new BigInteger("21039231"))); }
	//@Test public void huge1CC8ZeroRankStaticOptimalPruning() { assertThat(getZeroRank(makeHuge1CC8(), new StaticOptimalPruningOrderer()), is(new BigInteger("21039231"))); }
	//@Test public void huge1CC8ZeroRankDynamicHeuristicPruning() { assertThat(getZeroRank(makeHuge1CC8(), new DynamicHeuristicPruningOrderer()), is(new BigInteger("21039231"))); }
	// sequential:                  rank 21039231, finished in 14.29 s
	// dynamic opt pruning:         rank 21039231, finished in 1.35 m !!!! slow !!!!
	// static opt pruning:          rank 21039231, finished in 3.48 s
	// dynamic heuristic pruning:   rank 21039231, finished in 4.53 s

	private void checkEveryConf(Strand strand, ConfRanker.Orderer orderer) {

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
		ConfRanker ranker = new ConfRanker(confSpace, new RCs(confSpace), emat, orderer, (b) -> b.setTraditional());

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

	public BigInteger getZeroRank(Strand strand, ConfRanker.Orderer orderer) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		log("confs: %s", formatBig(new RCs(confSpace).getNumConformations()));

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
		ConfRanker ranker = new ConfRanker(confSpace, new RCs(confSpace), emat, orderer, (b) -> b.setTraditional());

		//ranker.reportProgress = true;

		Stopwatch stopwatch = new Stopwatch().start();
		BigInteger rank = ranker.getNumConfsAtMost(0.0);
		log("rank %s, finished in %s", rank, stopwatch.stop().getTime(2));

		return rank;
	}
}
