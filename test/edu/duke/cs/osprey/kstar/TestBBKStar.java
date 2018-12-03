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

import static edu.duke.cs.osprey.TestBase.TempFile;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.io.File;
import java.util.List;


public class TestBBKStar {

	public static class Results {
		public BBKStar bbkstar;
		public List<KStar.ScoredSequence> sequences;
	}

	public static Results runBBKStar(TestKStar.ConfSpaces confSpaces, int numSequences, double epsilon, String confdbPattern, int maxSimultaneousMutations,
									 boolean runMARKStar) {

		Parallelism parallelism = Parallelism.makeCpu(4);

		// how should we compute energies of molecules?
		try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.build()) {

			KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) ->
				String.format("TestBBKStar.assertSequence(results, \"%s\", %f, %f); // protein %s   ligand %s   complex %s",
					info.sequence.toString(Sequence.Renderer.ResType),
					info.kstarScore.lowerBoundLog10(),
					info.kstarScore.upperBoundLog10(),
					info.kstarScore.protein.toString(),
					info.kstarScore.ligand.toString(),
					info.kstarScore.complex.toString()
				);

			// configure BBK*
			KStar.Settings kstarSettings = new KStar.Settings.Builder()
				.setEpsilon(epsilon)
				.setStabilityThreshold(null)
				.setMaxSimultaneousMutations(maxSimultaneousMutations)
				.addScoreConsoleWriter(testFormatter)
				.build();
			BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
				.setNumBestSequences(numSequences)
				.setNumConfsPerBatch(8)
				.build();
			BBKStar bbkstar = new BBKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, kstarSettings, bbkstarSettings);
			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				// how should we define energies of conformations?
				info.confEcalcMinimized = new ConfEnergyCalculator.Builder(info.confSpace, ecalcMinimized)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, ecalcMinimized)
						.build()
						.calcReferenceEnergies()
					).build();

				// BBK* needs rigid energies too
				EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
					.setIsMinimizing(false)
					.build();
				ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(info.confEcalcMinimized, ecalcRigid);

				if(runMARKStar) {
					PartitionFunctionFactory pfuncFactory = new PartitionFunctionFactory(info.confSpace, info.confEcalcMinimized, info.id);
					pfuncFactory.setUseMARKStar(confEcalcRigid);
					info.pfuncFactory = pfuncFactory;
				}
				else {
					EnergyMatrix ematMinimized = new SimplerEnergyMatrixCalculator.Builder(info.confEcalcMinimized)
							.build()
							.calcEnergyMatrix();
					info.confSearchFactoryMinimized = (rcs) ->
							new ConfAStarTree.Builder(ematMinimized, rcs)
							.setTraditional()
							.build();
					EnergyMatrix ematRigid = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
							.build()
							.calcEnergyMatrix();
					info.confSearchFactoryRigid = (rcs) ->
							new ConfAStarTree.Builder(ematRigid, rcs)
							.setTraditional()
							.build();

				}

				// add the ConfDB file if needed
				if (confdbPattern != null) {
					info.confDBFile = new File(confdbPattern.replace("*", info.type.name().toLowerCase()));
				}
			}

			// run BBK*
			Results results = new Results();
			results.bbkstar = bbkstar;
			results.sequences = bbkstar.run();
			return results;
		}
	}

	public static Results runBBKStar(TestKStar.ConfSpaces confSpaces, int numSequences, double epsilon, String confdbPattern, int maxSimultaneousMutations) {
		return runBBKStar(confSpaces, numSequences, epsilon, confdbPattern, maxSimultaneousMutations, false);
	}

	@Test
	public void test2RL0() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.99;
		final int numSequences = 25;
		Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 10, true);

		assert2RL0(results, numSequences);
	}

	private void assert2RL0(Results results, int numSequences) {

		// K* bounds collected with e = 0.1 from original K* algo
		assertSequence(results, "PHE ASP GLU GLN PHE LYS ILE THR", 16.424659, 16.528464);
		assertSequence(results, "PHE ASP GLU ASN PHE LYS ILE THR", 16.114634, 16.216578);
		assertSequence(results, "PHE ASP GLU SER PHE LYS ILE THR", 15.970427, 16.078237);
		assertSequence(results, "PHE ASP GLU THR PHE LYS ILE SER", 15.881531, 15.985122);
		assertSequence(results, "PHE ASP GLU THR PHE LYS ILE THR", 15.315503, 15.394524);
		assertSequence(results, "PHE ASP GLU THR TYR LYS ILE THR", 15.029935, 15.108900);
		assertSequence(results, "PHE ASP GLU THR ILE LYS ILE THR", 15.005796, 15.081229);
		assertSequence(results, "PHE ASP GLU THR PHE LYS ILE ASN", 14.923698, 15.009355);
		assertSequence(results, "PHE ASP GLU THR VAL LYS ILE THR", 14.655568, 14.724370);
		assertSequence(results, "PHE ASP GLU THR LEU LYS ILE THR", 14.619748, 14.704462);
		assertSequence(results, "PHE ASP GLU THR PHE LYS VAL THR", 14.561234, 14.647761);
		assertSequence(results, "PHE ASP GLU THR PHE LYS LEU THR", 14.324171,14.405292);
		assertSequence(results, "PHE ASP GLU THR ALA LYS ILE THR", 14.159251, 14.225284);
		assertSequence(results, "PHE GLU GLU THR PHE LYS ILE THR", 14.056796, 14.148018);
		assertSequence(results, "PHE ASP GLU THR PHE LYS ALA THR", 13.987814, 14.064762);
		assertSequence(results, "PHE ASP ASP THR PHE LYS ILE THR", 13.412206, 13.489261);
		assertSequence(results, "ILE ASP GLU THR PHE LYS ILE THR", 12.844163, 12.936699);
		assertSequence(results, "VAL ASP GLU THR PHE LYS ILE THR", 12.612457, 12.677450);
		assertSequence(results, "LEU ASP GLU THR PHE LYS ILE THR", 12.336269, 12.417254);
		assertSequence(results, "ALA ASP GLU THR PHE LYS ILE THR", 11.778039, 11.840466);
		assertSequence(results, "TYR ASP GLU THR PHE LYS ILE THR", 11.550098, 11.633104);
		assertSequence(results, "PHE ASP GLU THR PHE ASP ILE THR", 10.805317, 10.871671);
		assertSequence(results, "PHE ASP GLU THR PHE GLU ILE THR", 10.012310, 10.079659);
		assertSequence(results, "PHE ASP GLU THR PHE LYS PHE THR", null, null);
		assertSequence(results, "PHE ASP GLU THR PHE LYS TYR THR", null, null);

		assertThat(results.sequences.size(), is(numSequences));
		assertDecreasingUpperBounds(results.sequences);
	}

	@Test
	public void test1GUA11() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make1GUA11();
		final double epsilon = 0.999999;
		final int numSequences = 6;
		Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1);

		// K* bounds collected with e = 0.1 from original K* algo
		assertSequence(results, "HIE ARG", 17.522258,17.636342);
		assertSequence(results, "HID VAL", 16.939674,17.014507);
		assertSequence(results, "HIE HIE", 16.833695,16.930972);
		assertSequence(results, "HIE HID", 16.659839,16.738627);
		assertSequence(results, "HIE LYS", 16.571112,16.683562);
		assertSequence(results, "HIE VAL", 16.474293,16.552681);

		assertThat(results.sequences.size(), is(numSequences));
		assertDecreasingUpperBounds(results.sequences);
	}

	@Test
	public void test2RL0WithConfDB() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.99;
		final int numSequences = 25;
		final String confdbPattern = "bbkstar.*.conf.db";

		try (TempFile proteinDBFile = new TempFile("bbkstar.protein.conf.db")) {
			try (TempFile ligandDBFile = new TempFile("bbkstar.ligand.conf.db")) {
				try (TempFile complexDBFile = new TempFile("bbkstar.complex.conf.db")) {

					// run with empty dbs
					Stopwatch sw = new Stopwatch().start();
					Results results = runBBKStar(confSpaces, numSequences, epsilon, confdbPattern, 1);
					assert2RL0(results, numSequences);
					System.out.println(sw.getTime(2));

					// the dbs should have stuff in them

					try (ConfDB confdb = new ConfDB(confSpaces.protein, proteinDBFile)) {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					try (ConfDB confdb = new ConfDB(confSpaces.ligand, ligandDBFile)) {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					try (ConfDB confdb = new ConfDB(confSpaces.complex, complexDBFile)) {
						assertThat(confdb.getNumSequences(), is((long)results.sequences.size()));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					}

					assertThat(proteinDBFile.exists(), is(true));
					assertThat(ligandDBFile.exists(), is(true));
					assertThat(complexDBFile.exists(), is(true));

					// run again with full dbs
					sw = new Stopwatch().start();
					Results results2 = runBBKStar(confSpaces, numSequences, epsilon, confdbPattern, 1);
					assert2RL0(results2, numSequences);
					System.out.println(sw.getTime(2));
				}
			}
		}
	}

	@Test
	public void only2RL0OneMutant() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0OnlyOneMutant();
		final double epsilon = 0.99;
		final int numSequences = 1;
		Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 1);

		assertThat(results.sequences.size(), is(1));
		assertThat(results.sequences.get(0).sequence.toString(Sequence.Renderer.AssignmentMutations), is("A193=VAL"));
	}

	@Test
	public void only2RL0SpaceWithoutWildType() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0SpaceWithoutWildType();
		final double epsilon = 0.99;
		final int numSequences = 2;
		Results results = runBBKStar(confSpaces, numSequences, epsilon, null, 2);

		assertThat(results.sequences.size(), is(2));
		assertThat(results.sequences.get(0).sequence.toString(Sequence.Renderer.AssignmentMutations), is("G654=thr A193=VAL"));
		assertThat(results.sequences.get(1).sequence.toString(Sequence.Renderer.AssignmentMutations), is("G654=VAL A193=VAL"));
	}

	private void assertSequence(Results results, String sequence, Double estKStarLowerLog10, Double estKStarUpperLog10) {

		// find the sequence
		for (KStar.ScoredSequence scoredSequence : results.sequences) {

			if (scoredSequence.sequence.toString(Sequence.Renderer.ResType).equals(sequence)) {

				// found it

				// check the K* bounds
				assertThat(scoredSequence.score, is(not(nullValue())));
				if (estKStarLowerLog10 != null && estKStarUpperLog10 != null) {

					// make sure these bounds contain the estimated bounds
					assertThat(scoredSequence.score.lowerBoundLog10(), lessThanOrEqualTo(estKStarLowerLog10));
					assertThat(scoredSequence.score.upperBoundLog10(), greaterThanOrEqualTo(estKStarUpperLog10));

				} else {
					assertThat(scoredSequence.score.score, is(nullValue()));
				}

				return;
			}
		}

		fail("sequence not found: " + sequence);
	}

	private void assertDecreasingUpperBounds(List<KStar.ScoredSequence> sequences) {

		double minUpperBoundLog10 = Double.POSITIVE_INFINITY;
		for (KStar.ScoredSequence sequence : sequences) {

			Double upperBoundLog10 = sequence.score.upperBoundLog10();
			if (upperBoundLog10 == null) {
				break;
			}
			assertThat(upperBoundLog10, lessThanOrEqualTo(minUpperBoundLog10));
			minUpperBoundLog10 = upperBoundLog10;
		}
	}
}
