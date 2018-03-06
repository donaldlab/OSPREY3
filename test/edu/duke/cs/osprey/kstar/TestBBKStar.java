package edu.duke.cs.osprey.kstar;

import static edu.duke.cs.osprey.TestBase.fileForWriting;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.Test;

import java.util.List;
import java.util.concurrent.atomic.AtomicReference;

public class TestBBKStar {

	public static class Results {
		public BBKStar bbkstar;
		public List<KStar.ScoredSequence> sequences;
	}

	public static Results runBBKStar(TestKStar.ConfSpaces confSpaces, int numSequences, double epsilon) {
		return runBBKStar(confSpaces, numSequences, epsilon, null);
	}

	public static Results runBBKStar(TestKStar.ConfSpaces confSpaces, int numSequences, double epsilon, String confdbPattern) {

		AtomicReference<Results> resultsRef = new AtomicReference<>(null);

		Parallelism parallelism = Parallelism.makeCpu(4);
		//Parallelism parallelism = Parallelism.make(4, 1, 8);

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.use((minimizingEcalc) -> {

				// how should we define energies of conformations?
				KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
					return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
						.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
							.build()
							.calcReferenceEnergies()
						).build();
				};

				// how should confs be ordered and searched?
				ConfSearchFactory confSearchFactory = (emat, rcs) -> {
					return new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();
				};

				// make a rigid energy calculator too
				EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
					.setIsMinimizing(false)
					.build();

				// run BBK*
				KStar.Settings kstarSettings = new KStar.Settings.Builder()
					.setEpsilon(epsilon)
					.setStabilityThreshold(null)
					.setMaxSimultaneousMutations(1)
					.addScoreConsoleWriter()
					.setConfDBPattern(confdbPattern)
					.build();
				BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
					.setNumBestSequences(numSequences)
					.setNumConfsPerBatch(8)
					.build();
				Results results = new Results();
				results.bbkstar = new BBKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, confSearchFactory, kstarSettings, bbkstarSettings);
				results.sequences = results.bbkstar.run();

				// pass back the ref
				resultsRef.set(results);
			});

		return resultsRef.get();
	}

	@Test
	public void test2RL0() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.99;
		final int numSequences = 25;
		Results results = runBBKStar(confSpaces, numSequences, epsilon);

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
		Results results = runBBKStar(confSpaces, numSequences, epsilon);

		// K* bounds collected with e = 0.1 from original K* algo
		assertSequence(results, "ILE ILE GLN HIE VAL TYR LYS ARG", 17.522258,17.636342);
		assertSequence(results, "ILE ILE GLN HID VAL TYR LYS VAL", 16.939674,17.014507);
		assertSequence(results, "ILE ILE GLN HIE VAL TYR LYS HIE", 16.833695,16.930972);
		assertSequence(results, "ILE ILE GLN HIE VAL TYR LYS HID", 16.659839,16.738627);
		assertSequence(results, "ILE ILE GLN HIE VAL TYR LYS LYS", 16.571112,16.683562);
		assertSequence(results, "ILE ILE GLN HIE VAL TYR LYS VAL", 16.474293,16.552681);

		assertThat(results.sequences.size(), is(numSequences));
		assertDecreasingUpperBounds(results.sequences);
	}

	@Test
	public void test2RL0WithConfDB() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.99;
		final int numSequences = 25;
		final String confdbPattern = "bbkstar.*.conf.db";

		fileForWriting("bbkstar.protein.conf.db", (proteinDBFile) -> {
			fileForWriting("bbkstar.ligand.conf.db", (ligandDBFile) -> {
				fileForWriting("bbkstar.complex.conf.db", (complexDBFile) -> {

					// run with empty dbs
					Stopwatch sw = new Stopwatch().start();
					Results results = runBBKStar(confSpaces, numSequences, epsilon, confdbPattern);
					assert2RL0(results, numSequences);
					System.out.println(sw.getTime(2));

					// the dbs should have stuff in them

					new ConfDB(confSpaces.protein, proteinDBFile).use((confdb) -> {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					new ConfDB(confSpaces.ligand, ligandDBFile).use((confdb) -> {
						assertThat(confdb.getNumSequences(), greaterThan(0L));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					new ConfDB(confSpaces.complex, complexDBFile).use((confdb) -> {
						assertThat(confdb.getNumSequences(), is((long)results.sequences.size()));
						for (Sequence sequence : confdb.getSequences()) {
							assertThat(confdb.getSequence(sequence).size(), greaterThan(0L));
						}
					});

					assertThat(proteinDBFile.exists(), is(true));
					assertThat(ligandDBFile.exists(), is(true));
					assertThat(complexDBFile.exists(), is(true));

					// run again with full dbs
					sw = new Stopwatch().start();
					Results results2 = runBBKStar(confSpaces, numSequences, epsilon, confdbPattern);
					assert2RL0(results2, numSequences);
					System.out.println(sw.getTime(2));
				});
			});
		});
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
