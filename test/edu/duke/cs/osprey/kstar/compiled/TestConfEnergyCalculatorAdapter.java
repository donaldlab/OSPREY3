package edu.duke.cs.osprey.kstar.compiled;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.Test;

import java.util.NoSuchElementException;
import java.util.function.Function;

import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;


public class TestConfEnergyCalculatorAdapter {

	@Test
	public void kstar() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"));

		KStarScoreWriter.Formatter testFormatter = (KStarScoreWriter.ScoreInfo info) -> {

			Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
				if (pfuncResult.status == PartitionFunction.Status.Estimated) {
					return String.format("%12e", pfuncResult.values.qstar.doubleValue());
				}
				return "null";
			};

			return String.format("assertSequence(result, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // protein %s ligand %s complex %s K* = %s",
				info.sequenceNumber,
				info.sequence.toString(Sequence.Renderer.ResType),
				formatPfunc.apply(info.kstarScore.protein),
				formatPfunc.apply(info.kstarScore.ligand),
				formatPfunc.apply(info.kstarScore.complex),
				info.kstarScore.protein,
				info.kstarScore.ligand,
				info.kstarScore.complex,
				info.kstarScore
			);
		};

		final double epsilon = 0.99;

		KStar.Settings settings = new KStar.Settings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(1)
			//.setShowPfuncProgress(true)
			.addScoreConsoleWriter(testFormatter)
			.build();
		TestKStar.Result result = new TestKStar.Result();
		result.kstar = new KStar(chainG, chainA, complex, settings);

		try (ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(4);

			for (KStar.ConfSpaceInfo info : result.kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
				boolean minimize = true;
				ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace, tasks);

				SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
					.setMinimize(minimize)
					.build()
					.calc();

				EnergyMatrix emat = new EmatCalculator.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(minimize)
					.build()
					.calc();

				info.confEcalc = new ConfEnergyCalculatorAdapter.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(minimize)
					.build();

				Function<RCs,ConfSearch> confSearchFactory = (rcs) ->
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();

				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalc,
					confSearchFactory.apply(rcs),
					confSearchFactory.apply(rcs),
					rcs.getNumConformations()
				);
			}

			// run K*
			result.scores = result.kstar.run(tasks);

			// cleanup
			for (KStar.ConfSpaceInfo info : result.kstar.confSpaceInfos()) {
				((ConfEnergyCalculatorAdapter)info.confEcalc).confEcalc.close();
			}
		}

		// make sure the K* bounds are correct
		assertKStar(result, "PHE LYS ILE SER PHE ASP GLU", 19.623811, 19.748049);
		assertKStar(result, "PHE LYS ILE THR PHE ASP GLU", 19.203386, 19.323441);
		assertKStar(result, "TYR LYS ILE THR PHE ASP GLU", 19.228177, 19.345113);
		assertKStar(result, "ILE LYS ILE THR PHE ASP GLU", 18.982370, 19.103969);
		assertKStar(result, "PHE LYS ILE THR PHE GLU GLU", 18.902824, 19.023228);
		assertKStar(result, "PHE LYS ILE ASN PHE ASP GLU", 18.802481, 18.910993);
		assertKStar(result, "PHE LYS VAL THR PHE ASP GLU", 18.740128, 18.862470);
		assertKStar(result, "VAL LYS ILE THR PHE ASP GLU", 18.595759, 18.710352);
		assertKStar(result, "LEU LYS ILE THR PHE ASP GLU", 18.202775, 18.320629);
		assertKStar(result, "PHE LYS ALA THR PHE ASP GLU", 18.119919, 18.241689);
		assertKStar(result, "ALA LYS ILE THR PHE ASP GLU", 18.087688, 18.198495);
		assertKStar(result, "PHE LYS ILE THR LEU ASP GLU", 17.514978, 17.632124);
		assertKStar(result, "PHE LYS ILE THR PHE ASP ASP", 17.009586, 17.121840);
		assertKStar(result, "PHE LYS ILE THR ILE ASP GLU", 16.961519, 17.078603);
		assertKStar(result, "PHE LYS LEU THR PHE ASP GLU", 16.918166, 17.034486);
		assertKStar(result, "PHE LYS ILE THR VAL ASP GLU", 16.459856, 16.565505);
		assertKStar(result, "PHE LYS ILE THR TYR ASP GLU", 15.884875, 16.004093);
		assertKStar(result, "PHE LYS ILE THR ALA ASP GLU", 15.567365, 15.671388);
		assertKStar(result, "PHE ASP ILE THR PHE ASP GLU", 14.300678, 14.411186);
		assertKStar(result, "PHE GLU ILE THR PHE ASP GLU", 13.555692, 13.652010);
		assertKStar(result, "PHE LYS PHE THR PHE ASP GLU", Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
		assertKStar(result, "PHE LYS TYR THR PHE ASP GLU", Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);

		/* NOTE:
			These results were computed with e=0.1

			At the time of test writing, the ranking of the sequences computed by the new code matched the ranking
			produced by the existing K* code, but the K* scores computed by this code are quite different.
			That's actually ok, as long as the rankings match.

			Ideally, we'd test the rankings themselves instead of the K* scores here, but
			to accurately rank the sequences requires computing with a small epsilon (eg < 0.1).
			Since computing with small epsilons takes forever, that don't make a very good test case.

			So the next best thing we can do is compute with a large epsilon, and make sure the K* scores
			still match. If the K* scores still match, that means the rankings must match too and we're good.

			Even with the large epsilon, this test still takes about a minute on my laptop.
		*/
	}

	private void assertKStar(TestKStar.Result result, String seqstr, double minKStar, double maxKStar) {

		MathTools.DoubleBounds exp = new MathTools.DoubleBounds(minKStar, maxKStar);

		// find the sequence in the result
		KStar.ScoredSequence seq = result.scores.stream()
			.filter(scoredSeq -> scoredSeq.sequence.toString(Sequence.Renderer.ResType).equals(seqstr))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("didn't find seuqnece: " + seqstr));

		MathTools.DoubleBounds obs = new MathTools.DoubleBounds(
			seq.score.lowerBoundLog10(),
			seq.score.upperBoundLog10()
		);

		// the two K* bounds should overlap
		assertThat(
			String.format("exp: %s\nobs: %s\n", exp, obs),
			obs.intersects(exp), is(true)
		);
	}
}
