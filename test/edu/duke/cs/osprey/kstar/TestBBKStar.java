package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicReference;

public class TestBBKStar {

	public static class Results {
		public BBKStar bbkstar;
		public BBKStar.ScoredSequence wildtype;
		public List<BBKStar.ScoredSequence> mutants;
	}

	public static Results runBBKStar(TestKStar.ConfSpaces confSpaces, double epsilon) {

		AtomicReference<Results> resultsRef = new AtomicReference<>(null);

		Parallelism parallelism = Parallelism.makeCpu(4);

		// how should we compute energies of molecules?
		new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.use((ecalc) -> {

				// how should we define energies of conformations?
				KStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
					return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
						.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
							.build()
							.calcReferenceEnergies()
						).build();
				};

				// how should confs be ordered and searched?
				ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						.setTraditional()
						.build();
				};

				// run BBK*
				BBKStar.Settings settings = new BBKStar.Settings.Builder()
					.setEpsilon(epsilon)
					.setMaxSimultaneousMutations(1)
					.setNumBestSequences(24)
					//.setShowPfuncProgress(true) // TEMP
					.setNumConfsPerBatch(8)
					.build();
				Results results = new Results();
				results.bbkstar = new BBKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, ecalc, confEcalcFactory, confSearchFactory, settings);
				results.wildtype = results.bbkstar.calcWildType();
				results.mutants = results.bbkstar.calcMutants();

				// pass back the ref
				resultsRef.set(results);
			});

		return resultsRef.get();
	}

	// NOTE: this test takes ~7 minutes to finish on GPU hardware, so it's disabled by default
	// TODO: need a shorter test for routine regression testing
	// sequence order won't necessarily be preserved with a larger epsilon, so can't save time
	// TEMP: test enabled
	@Test
	public void test2RL0() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.1;
		Results results = runBBKStar(confSpaces, epsilon);

		assertThat(results.wildtype, is(not(nullValue())));
		assertThat(results.mutants.size(), is(24));

		// scored collected with e = 0.1 from original K* algo
		assertSequence(results.wildtype,       "PHE ASP GLU THR PHE LYS ILE THR", 15.340604, epsilon);
		assertSequence(results.mutants.get(0), "PHE ASP GLU GLN PHE LYS ILE THR", 16.476261, epsilon);
		assertSequence(results.mutants.get(1), "PHE ASP GLU ASN PHE LYS ILE THR", 16.163977, epsilon);
		assertSequence(results.mutants.get(2), "PHE ASP GLU SER PHE LYS ILE THR", 16.026221, epsilon);
		assertSequence(results.mutants.get(3), "PHE ASP GLU THR PHE LYS ILE SER", 15.933782, epsilon);
		assertSequence(results.mutants.get(4), "PHE ASP GLU THR TYR LYS ILE THR", 15.057543, epsilon);
		assertSequence(results.mutants.get(5), "PHE ASP GLU THR ILE LYS ILE THR", 15.029971, epsilon);
		assertSequence(results.mutants.get(6), "PHE ASP GLU THR PHE LYS ILE ASN", 14.953614, epsilon);
		assertSequence(results.mutants.get(7), "PHE ASP GLU THR VAL LYS ILE THR", 14.673854, epsilon);
		assertSequence(results.mutants.get(8), "PHE ASP GLU THR LEU LYS ILE THR", 14.648011, epsilon);
		assertSequence(results.mutants.get(9), "PHE ASP GLU THR PHE LYS VAL THR", 14.592144, epsilon);
		assertSequence(results.mutants.get(10), "PHE ASP GLU THR ALA LYS ILE THR", 14.175944, epsilon);
		assertSequence(results.mutants.get(11), "PHE GLU GLU THR PHE LYS ILE THR", 14.094045, epsilon);
		assertSequence(results.mutants.get(12), "PHE ASP GLU THR PHE LYS ALA THR", 14.008732, epsilon);
		assertSequence(results.mutants.get(13), "PHE ASP ASP THR PHE LYS ILE THR", 13.438073, epsilon);
		assertSequence(results.mutants.get(14), "ILE ASP GLU THR PHE LYS ILE THR", 12.882685, epsilon);
		assertSequence(results.mutants.get(15), "VAL ASP GLU THR PHE LYS ILE THR", 12.628520, epsilon);
		assertSequence(results.mutants.get(16), "LEU ASP GLU THR PHE LYS ILE THR", 12.366343, epsilon);
		assertSequence(results.mutants.get(17), "ALA ASP GLU THR PHE LYS ILE THR", 11.792450, epsilon);
		assertSequence(results.mutants.get(18), "TYR ASP GLU THR PHE LYS ILE THR", 11.580947, epsilon);
		assertSequence(results.mutants.get(19), "PHE ASP GLU THR PHE ASP ILE THR", 10.823063, epsilon);
		assertSequence(results.mutants.get(20), "PHE ASP GLU THR PHE GLU ILE THR", 10.031569, epsilon);

		// the unscorable sequences could happen in any order
		Set<String> unscorableSequences = new HashSet<>();
		unscorableSequences.add(results.mutants.get(21).sequence.toString());
		unscorableSequences.add(results.mutants.get(22).sequence.toString());
		unscorableSequences.add(results.mutants.get(23).sequence.toString());
		assertThat(unscorableSequences, containsInAnyOrder(
			"PHE ASP GLU THR PHE LYS LEU THR",
			"PHE ASP GLU THR PHE LYS PHE THR",
			"PHE ASP GLU THR PHE LYS TYR THR"
		));
	}

	private void assertSequence(BBKStar.ScoredSequence scoredSequence, String sequence, Double estKstarLog10, double pfuncEpsilon) {

		// check the sequence
		assertThat(scoredSequence.sequence.toString(), is(sequence));

		// check the K* score
		Double kstarScoreLog10 = scoredSequence.score.scoreLog10();
		if (estKstarLog10 == null) {
			assertThat(kstarScoreLog10, is(nullValue()));
		} else {
			assertThat(kstarScoreLog10, is(not(nullValue())));

			// combine the pfunc epsilons to get an epsilon value effective for K* scores
			// (and convert to log10 space)
			double kstarEpsilonLog10 = Math.log10(1.0/(1.0 - pfuncEpsilon));

			// we don't know the true K*, but we can bound it based on the estimate given
			double minKStarLog10 = estKstarLog10;
			double maxKStarLog10 = minKStarLog10 + 2*kstarEpsilonLog10;

			// make sure the observed score is within epsilon of the bound
			assertThat(kstarScoreLog10, greaterThanOrEqualTo(minKStarLog10 - kstarEpsilonLog10));
			assertThat(kstarScoreLog10, lessThanOrEqualTo(maxKStarLog10 + kstarEpsilonLog10));
		}
	}
}
