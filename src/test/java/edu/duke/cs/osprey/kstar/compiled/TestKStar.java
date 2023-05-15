package edu.duke.cs.osprey.kstar.compiled;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.ematrix.compiled.ErefCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculatorAdapter;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.FileTools;
import edu.duke.cs.osprey.tools.Log;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import org.junit.jupiter.api.Test;

import java.time.Duration;
import java.util.List;
import java.util.function.Function;

import static org.hamcrest.Matchers.is;
import static org.hamcrest.MatcherAssert.assertThat;


public class TestKStar {

	private boolean generate = false;
	private double epsilon = Double.NaN;

	public static class Result {
		public KStar kstar;
		public List<KStar.ScoredSequence> scores;
	}


	private void test2RL0(PosInterDist posInterDist) {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"));

		final double epsilon = 0.99;
		var result = run(chainG, chainA, complex, posInterDist, epsilon);

		// check the results
		// answers computed with epsilon=0.1, on 48 threads, with a 5 minute timeout for each pfunc
		assertSequence(result,   0, "PHE LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "4.697497e+934", "5.044775e+934", "1.695936e+1117", "1.876755e+1117");
		assertSequence(result,   1, "PHE LYS ILE THR PHE ASP ASP", "3.792600e+117", "4.129520e+117", "4.697497e+934", "5.044775e+934", "2.404320e+1113", "2.596913e+1113");
		assertSequence(result,   2, "PHE LYS ILE THR PHE GLU GLU", "2.201263e+119", "2.245939e+119", "4.697497e+934", "5.044775e+934", "1.015415e+1117", "1.115463e+1117");
		assertSequence(result,   3, "PHE LYS ILE THR ALA ASP GLU", "4.565707e+116", "4.565707e+116", "4.697497e+934", "5.044775e+934", "1.040126e+1111", "1.111272e+1111");
		assertSequence(result,   4, "PHE LYS ILE THR ILE ASP GLU", "3.062507e+116", "3.070559e+116", "4.697497e+934", "5.044775e+934", "1.626136e+1112", "1.797470e+1112");
		assertSequence(result,   5, "PHE LYS ILE THR LEU ASP GLU", "4.078772e+117", "4.123589e+117", "4.697497e+934", "5.044775e+934", "8.573269e+1113", "9.389148e+1113");
		assertSequence(result,   6, "PHE LYS ILE THR TYR ASP GLU", "2.166241e+119", "2.179750e+119", "4.697497e+934", "5.044775e+934", "1.079470e+1114", "1.197624e+1114");
		assertSequence(result,   7, "PHE LYS ILE THR VAL ASP GLU", "7.693108e+115", "7.804112e+115", "4.697497e+934", "5.044775e+934", "1.383324e+1111", "1.497185e+1111");
		assertSequence(result,   8, "PHE LYS ILE ASN PHE ASP GLU", "1.599167e+119", "1.603392e+119", "4.120907e+933", "4.131389e+933", "4.525219e+1115", "4.653164e+1115");
		assertSequence(result,   9, "PHE LYS ILE SER PHE ASP GLU", "1.599167e+119", "1.603392e+119", "2.009793e+934", "2.215483e+934", "1.720875e+1117", "1.863201e+1117");
		assertSequence(result,  10, "PHE LYS ALA THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "4.112727e+931", "4.392889e+931", "1.326162e+1113", "1.473341e+1113");
		assertSequence(result,  11, "PHE LYS LEU THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "4.639921e+903", "4.946852e+903", "8.850462e+1083", "9.648729e+1083");
		assertSequence(result,  12, "PHE LYS PHE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "1.250916e+930", "1.343667e+930", "1.214005e+885", "2.684640e+905");
		assertSequence(result,  13, "PHE LYS TYR THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "2.225535e+930", "2.390111e+930", "2.176774e+828", "1.104128e+864");
		assertSequence(result,  14, "PHE LYS VAL THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "8.928763e+932", "9.483391e+932", "9.999413e+1114", "1.103053e+1115");
		assertSequence(result,  15, "PHE ASP ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "2.916204e+926", "2.949550e+926", "1.149242e+1104", "1.251620e+1104");
		assertSequence(result,  16, "PHE GLU ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "1.102772e+927", "1.104368e+927", "7.773140e+1103", "8.353264e+1103");
		assertSequence(result,  17, "ALA LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "1.676580e+931", "1.701650e+931", "4.094726e+1112", "4.535013e+1112");
		assertSequence(result,  18, "ILE LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "9.521737e+932", "1.056179e+933", "1.999694e+1115", "2.209177e+1115");
		assertSequence(result,  19, "LEU LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "7.164489e+932", "7.452615e+932", "1.942237e+1114", "2.148226e+1114");
		assertSequence(result,  20, "TYR LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "7.725090e+934", "8.475823e+934", "2.815613e+1117", "3.103106e+1117");
		assertSequence(result,  21, "VAL LYS ILE THR PHE ASP GLU", "1.599167e+119", "1.603392e+119", "2.647316e+932", "2.727974e+932", "2.073171e+1114", "2.296308e+1114");
	}
	@Test public void test2RL0_DesmetEtAl1992() { test2RL0(PosInterDist.DesmetEtAl1992); }
	// tragically, calculating the energy matrices for the tighter bounds takes longer than the time limit here, so this doesn't make a good test case
	/* @Test */ public void test2RL0_TighterBounds() { test2RL0(PosInterDist.TighterBounds); }

	// TODO: add more test cases here?

	private Result run(ConfSpace protein, ConfSpace ligand, ConfSpace complex, PosInterDist posInterDist, double epsilon) {

		var settings = new KStar.Settings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(1)
			.setPfuncTimeout(Duration.ofSeconds(10));

		if (generate) {

			// override the settings to generate the expected vaules
			settings.setEpsilon(this.epsilon);
			settings.setPfuncTimeout(Duration.ofMinutes(5));

			settings.addScoreConsoleWriter(info -> {

				Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
					var bounds = pfuncResult.values.calcBounds();
					return String.format("\"%-12s\", \"%-12s\"",
						Log.formatBigEngineering(bounds.lower),
						Log.formatBigEngineering(bounds.upper)
					);
				};

				return String.format("assertSequence(result, %3d, \"%s\", %s, %s, %s);",
					info.sequenceNumber,
					info.sequence.toString(Sequence.Renderer.ResType),
					formatPfunc.apply(info.kstarScore.protein),
					formatPfunc.apply(info.kstarScore.ligand),
					formatPfunc.apply(info.kstarScore.complex)
				);
			});

		} else {

			settings.addScoreConsoleWriter(info ->
				String.format("finished sequence %d/%d [%s]",
					info.sequenceNumber + 1,
					info.numSequences,
					info.sequence.toString(Sequence.Renderer.ResType)
				)
			);
		}

		KStar kstar = new KStar(protein, ligand, complex, settings.build());

		var parallelism = Parallelism.makeCpu(4);
		try (var tasks = parallelism.makeTaskExecutor()) {

			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				ConfEnergyCalculator ecalc = ConfEnergyCalculator.makeBest(confSpace);

				SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
					.build()
					.calc(tasks);

				EnergyMatrix emat = new EmatCalculator.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.build()
					.calc(tasks);

				info.confEcalc = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.build();

				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalc,
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build(),
					rcs.getNumConformations()
				).setPreciseBcalc(true);
			}

			// run K*
			var result = new Result();
			result.kstar = kstar;
			result.scores = kstar.run(tasks);
			return result;

		} finally {

			// cleanup
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
				if (info.confEcalc != null) {
					((ConfEnergyCalculatorAdapter)info.confEcalc).confEcalc.close();
				}
			}
		}
	}

	private void assertSequence(Result result, int index, String seqstr, String proteinMin, String proteinMax, String ligandMin, String ligandMax, String complexMin, String complexMax) {

		var proteinExp = new BigDecimalBounds(proteinMin, proteinMax);
		var ligandExp = new BigDecimalBounds(ligandMin, ligandMax);
		var complexExp = new BigDecimalBounds(complexMin, complexMax);

		KStar.ScoredSequence seq = result.scores.get(index);
		assertThat(seq.sequence.toString(Sequence.Renderer.ResType), is(seqstr));

		var proteinObs = seq.score.protein.values.calcBounds();
		var ligandObs = seq.score.ligand.values.calcBounds();
		var complexObs = seq.score.complex.values.calcBounds();

		// the two pfunc bounds should overlap for each state
		assertThat(
			String.format("\n\texp=%s\n\tobs=%s", proteinExp, proteinObs),
			proteinObs.intersects(proteinExp),
			is(true)
		);
		assertThat(
			String.format("\n\texp=%s\n\tobs=%s", ligandExp, ligandObs),
			ligandObs.intersects(ligandExp),
			is(true)
		);
		assertThat(
			String.format("\n\texp=%s\n\tobs=%s", complexExp, complexObs),
			complexObs.intersects(complexExp),
			is(true)
		);
	}

	public static void main(String[] args) {

		var test = new TestKStar();
		test.generate = true;

		test.epsilon = 0.1;
		test.test2RL0_TighterBounds();
	}
}
