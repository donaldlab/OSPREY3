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

import java.math.BigInteger;
import java.util.Arrays;


public class TestConfRanker {

	public static Strand makeTinyDiscrete1CC8() {

		// 8 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A7").setLibraryRotamers("GLY", "ALA");
		strand.flexibility.get("A9").setLibraryRotamers("GLY", "ALA");
		return strand;
	}
	@Test public void tinyDiscrete1CC8() { checkEveryConf(makeTinyDiscrete1CC8()); }

	public static Strand makeSmall1CC8() {

		// 2268 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A5").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		strand.flexibility.get("A7").setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		return strand;
	}
	@Test public void small1CC8() { checkEveryConf(makeSmall1CC8()); }

	public static Strand makeMedium1CC8() {

		// 3.62e6 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	@Test public void medium1CC8() { assertThat(getZeroRank(makeMedium1CC8()), is(new BigInteger("40306"))); }

	public static Strand makeLarge1CC8() {

		// 3.62e6 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	@Test public void large1CC8() { assertThat(getZeroRank(makeLarge1CC8()), is(new BigInteger("1034629"))); }

	public static Strand makeHuge1CC8() {

		// 3.86e9 confs
		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14")) {
			strand.flexibility.get(resNum).setLibraryRotamers(Strand.WildType).addWildTypeRotamers().setContinuous();
		}
		return strand;
	}
	@Test public void huge1CC8() { assertThat(getZeroRank(makeHuge1CC8()), is(new BigInteger("21039231"))); }

	// benchmark data for various position ordering heuristics attempted in the past:
	// sequential:                    rank 21039231, finished in 14.29 s
	// dynamic opt pruning:           rank 21039231, finished in 1.35 m !!!! slow !!!!
	// static opt pruning:            rank 21039231, finished in 3.48 s
	// dynamic heuristic:             rank 21039231, finished in 4.53 s
	// optimized dynamic heuristic:   rank 21039231, finished in 1.16 s

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
		ConfRanker ranker = new ConfRanker.Builder(confSpace, emat)
			.build();

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

	public BigInteger getZeroRank(Strand strand) {

		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		//log("confs: %s", formatBig(new RCs(confSpace).getNumConformations()));

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
		ConfRanker ranker = new ConfRanker.Builder(confSpace, emat)
			.build();

		/*
		ranker.reportProgress = true;
		Stopwatch stopwatch = new Stopwatch().start();
		BigInteger rank = ranker.getNumConfsAtMost(0.0);
		log("rank %s, finished in %s", rank, stopwatch.stop().getTime(2));
		return rank;
		*/

		return ranker.getNumConfsAtMost(0.0);
	}
}
