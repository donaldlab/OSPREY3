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

package edu.duke.cs.osprey.pruning;

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;

public class TestSimpleDEE extends TestBase {

	/**
	 * use this to generate encoded pmats for unit tests
	 * (you might need to hack up PruningControl and Pruner though to isolate different kinds of pruning)
	 */
	public static void main(String[] args) {

		initDefaultEnvironment();

		// get an old-school search problem so we can use the original DEE code
		TestBase.ResidueFlexibility resFlex = new TestBase.ResidueFlexibility();
		resFlex.addFlexible("3 4 6");
		resFlex.sortPositions();
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

		// calc an old-school energy matrix
		EnergyMatrixCalculator emcalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, search.useERef, search.addResEntropy);
		emcalc.calcPEM();
		search.emat = emcalc.getEMatrix();
		System.out.println(search.emat.toString());

		// do old-school DEE
		double pruningInterval = 100.0; // I0 + Ew
		boolean typeDep = false; // default = false
		double boundsThreshold = 0; // default = 100, but never used in DEE as far as I can tell
		int algOption = 3; // default = 1
		boolean useFlags = true; // default = true
		boolean useTriples = false; // default = false
		boolean preDACS = false; // default = false
		double stericThreshold = Double.POSITIVE_INFINITY; // default = 100
		search.pruneMat = new PruningMatrix(search.confSpace, 0);
		new PruningControl(search,
			pruningInterval, typeDep, boundsThreshold, algOption, useFlags,
			useTriples, preDACS, useEpic, useTupleExpansion, stericThreshold
		).prune();

		// update pruned pairs with singles info
		search.pruneMat.prunePairsFromSingles();

		System.out.println(search.pruneMat);
		System.out.println(writePmat(search.pruneMat));
	}

	public static String writePmat(PruningMatrix pmat) {

		StringBuilder buf = new StringBuilder();

		Consumer<Boolean> writer = (val) -> {
			buf.append(val ? '+' : '-');
		};

		for (int res1=0; res1<pmat.getNumPos(); res1++) {
			int m1 = pmat.getNumConfAtPos(res1);
			for (int i1=0; i1<m1; i1++) {
				writer.accept(pmat.getOneBody(res1, i1));
				for (int res2=0; res2<res1; res2++) {
					int m2 = pmat.getNumConfAtPos(res2);
					for (int i2=0; i2<m2; i2++) {
						writer.accept(pmat.getPairwise(res1, i1, res2, i2));
					}
				}
			}
		}

		return buf.toString();
	}

	public static PruningMatrix readPmat(SimpleConfSpace confSpace, String encoded) {
		PruningMatrix pmat = new PruningMatrix(confSpace);
		pmat.fill(new Iterator<Boolean>() {

			int pos = 0;

			@Override
			public boolean hasNext() {
				return pos < encoded.length();
			}

			@Override
			public Boolean next() {
				boolean val = encoded.charAt(pos) == '+';
				pos++;
				return val;
			}
		});

		// SimpleDEE prunes pairs if either of the corresponding singles are pruned
		// so do the same for the old pmats, to be comparable
		pmat.prunePairsFromSingles();

		return pmat;
	}

	private PruningMatrix calcPmat(SimpleConfSpace confSpace, Consumer<SimpleDEE.Runner> func) {

		AtomicReference<PruningMatrix> pmat = new AtomicReference<>(null);

		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.use((ecalc) -> {

				// calc an energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
					.build()
					.calcEnergyMatrix();

				// calc a pruning matrix with DEE
				SimpleDEE.Runner runner = new SimpleDEE.Runner();
				func.accept(runner);
				runner.setShowProgress(true);
				pmat.set(runner.run(confSpace, emat));
			});

		return pmat.get();
	}

	public static SimpleConfSpace make1CC8_3Pos() {
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		strand.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A6").setLibraryRotamers(Strand.WildType);
		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	@Test
	public void test1CC8_3Pos_SinglesThreshold() {
		SimpleConfSpace confSpace = make1CC8_3Pos();
		PruningMatrix pmat = calcPmat(confSpace, (runner) -> {
			runner.setSinglesThreshold(100.0);
			runner.setPairsThreshold(null);
		});
		assertThat(pmat, is(readPmat(confSpace, "----+------------------------------------------------------------------++++++++++++++++++++++++++++++++----------------++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")));
	}

	@Test
	public void test1CC8_3Pos_PairsThreshold() {
		SimpleConfSpace confSpace = make1CC8_3Pos();
		PruningMatrix pmat = calcPmat(confSpace, (runner) -> {
			runner.setSinglesThreshold(null);
			runner.setPairsThreshold(0.0);
		});
		assertThat(pmat, is(readPmat(confSpace, "-------------+--------+--------+--------+--------+--------+--------+----++++++++------------+-----------+-+++-+---------+-+++-+---------+-+++-++------------+--------------++-----------+-+++-++-------")));
	}

	@Test
	public void test1CC8_3Pos_SinglesGoldstein() {
		SimpleConfSpace confSpace = make1CC8_3Pos();
		PruningMatrix pmat = calcPmat(confSpace, (runner) -> {
			runner.setSinglesThreshold(null);
			runner.setPairsThreshold(null);
			runner.setSinglesGoldsteinDiffThreshold(100.0);
		});
		assertThat(pmat, is(readPmat(confSpace, "----+---+++++++++------------------------------------------------------++++++++++++++++++++++++++++++++----------------++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")));
	}

	@Test
	public void test1CC8_3Pos_PairsGoldstein() {
		SimpleConfSpace confSpace = make1CC8_3Pos();
		PruningMatrix pmat = calcPmat(confSpace, (runner) -> {
			runner.setSinglesThreshold(null);
			runner.setPairsThreshold(null);
			runner.setPairsGoldsteinDiffThreshold(100.0);
		});
		assertThat(pmat, is(readPmat(confSpace, "---------++++++++--+--+-----+--+-----+--+--------+--------+--------+----+++++++++++++++-+++++++++++++++-----+---+-------+++++++++++++++-+++++++++++++++-+++++++++++++++-+++++++++++++++-+++++++++++++++")));
	}

	@Test
	public void prunedConfsOutsideInterval() {

		// check that pruned confs are actually greater than the min-bound conf plus the interval
		final double interval = 4.0;

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		for (String resNum : Arrays.asList("A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10")) { // keep small enough we can exhuastively enumerate conf space
			strand.flexibility.get(resNum).setLibraryRotamers("ALA", "VAL", "LEU", "ILE");
		}
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		EnergyMatrix emat;
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			emat = new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
				.build()
				.calcEnergyMatrix();
		}

		PruningMatrix pmat = new SimpleDEE.Runner()
			.setThreshold(null)
			.setSinglesGoldsteinDiffThreshold(interval)
			//.setPairsGoldsteinDiffThreshold(interval)
			//.setTriplesGoldsteinDiffThreshold(interval)
			.setShowProgress(true)
			.run(confSpace, emat);

		ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
			.setTraditional()
			.build();

		ConfSearch.ScoredConf minScoreConf = astar.nextConf();

		// TEMP
		log("min score conf: %.4f", minScoreConf.getScore());

		while (true) {
			ConfSearch.ScoredConf conf = astar.nextConf();
			if (conf == null) {
				break;
			}

			if (pmat.isPruned(new RCTuple(conf.getAssignments()))) {
				double diff = conf.getScore() - minScoreConf.getScore();
				log("min score pruned conf: %.4f   (diff: %.4f)", conf.getScore(), diff);
				assertThat(diff, greaterThan(interval));
				break;
			}
		}
	}
}
