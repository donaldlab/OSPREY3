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
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.Timeout;

import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;


public class TestKStar {

	// buggy versions of these tests sometimes take a very long time to run
	@Rule
	public Timeout globalTimeout = new Timeout(5, TimeUnit.MINUTES);

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
		// answers computed with epsilon=0.1 using the GPU ecalc
		assertSequence(result,   0, "PHE LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "4.723506e+934", "5.144046e+934", "3.712483e+1106", "4.124440e+1106");
		assertSequence(result,   1, "PHE LYS ILE THR PHE ASP ASP", "3.698526e+117", "3.964339e+117", "4.723506e+934", "5.144046e+934", "8.815333e+1104", "9.792659e+1104");
		assertSequence(result,   2, "PHE LYS ILE THR PHE GLU GLU", "2.183036e+119", "2.374596e+119", "4.723506e+934", "5.144046e+934", "5.219978e+1106", "5.799442e+1106");
		assertSequence(result,   3, "PHE LYS ILE THR ALA ASP GLU", "4.452034e+116", "4.715925e+116", "4.723506e+934", "5.144046e+934", "1.053552e+1104", "1.170024e+1104");
		assertSequence(result,   4, "PHE LYS ILE THR ILE ASP GLU", "2.993563e+116", "3.234977e+116", "4.723506e+934", "5.144046e+934", "7.158340e+1103", "7.952307e+1103");
		assertSequence(result,   5, "PHE LYS ILE THR LEU ASP GLU", "4.011891e+117", "4.363522e+117", "4.723506e+934", "5.144046e+934", "9.591381e+1104", "1.065557e+1105");
		assertSequence(result,   6, "PHE LYS ILE THR TYR ASP GLU", "2.109162e+119", "2.287207e+119", "4.723506e+934", "5.144046e+934", "5.048081e+1106", "5.608080e+1106");
		assertSequence(result,   7, "PHE LYS ILE THR VAL ASP GLU", "7.490372e+115", "7.946553e+115", "4.723506e+934", "5.144046e+934", "1.782366e+1103", "1.979794e+1103");
		assertSequence(result,   8, "PHE LYS ILE ASN PHE ASP GLU", "1.545472e+119", "1.684915e+119", "3.918142e+933", "4.205637e+933", "1.533457e+1105", "1.702591e+1105");
		assertSequence(result,   9, "PHE LYS ILE SER PHE ASP GLU", "1.545472e+119", "1.684915e+119", "2.063374e+934", "2.281592e+934", "2.232865e+1106", "2.480699e+1106");
		assertSequence(result,  10, "PHE LYS ALA THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "3.980326e+931", "4.382659e+931", "1.601785e+1103", "1.779294e+1103");
		assertSequence(result,  11, "PHE LYS LEU THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "4.687707e+903", "5.131922e+903", "3.218183e+1075", "3.574845e+1075");
		assertSequence(result,  12, "PHE LYS PHE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "1.274839e+930", "1.388108e+930", "1.339810e+1103", "1.488344e+1103");
		assertSequence(result,  13, "PHE LYS TYR THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "2.274307e+930", "2.497399e+930", "3.926699e+1103", "4.362458e+1103");
		assertSequence(result,  14, "PHE LYS VAL THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "8.756688e+932", "9.644181e+932", "5.116884e+1104", "5.684792e+1104");
		assertSequence(result,  15, "PHE ASP ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "2.756439e+926", "2.983265e+926", "7.913709e+1095", "8.787334e+1095");
		assertSequence(result,  16, "PHE GLU ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "1.070749e+927", "1.114131e+927", "1.578784e+1096", "1.752743e+1096");
		assertSequence(result,  17, "ALA LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "1.622001e+931", "1.747574e+931", "9.186497e+1101", "1.020473e+1102");
		assertSequence(result,  18, "ILE LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "9.889412e+932", "1.086031e+933", "4.469246e+1104", "4.964611e+1104");
		assertSequence(result,  19, "LEU LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "6.964435e+932", "7.652128e+932", "4.328846e+1103", "4.808722e+1103");
		assertSequence(result,  20, "TYR LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "7.924623e+934", "8.691902e+934", "6.128211e+1106", "6.807914e+1106");
		assertSequence(result,  21, "VAL LYS ILE THR PHE ASP GLU", "1.545472e+119", "1.684915e+119", "2.563908e+932", "2.780140e+932", "4.652633e+1103", "5.168006e+1103");
	}
	@Test public void test2RL0_DesmetEtAl1992() { test2RL0(PosInterDist.DesmetEtAl1992); }
	// tragically, calculating the energy matrices for the tighter bounds takes longer than the time limit here, so this doesn't make a good test case
	/* @Test */ public void test2RL0_TighterBounds() { test2RL0(PosInterDist.TighterBounds); }

	// TODO: add more test cases here?


	private Result run(ConfSpace protein, ConfSpace ligand, ConfSpace complex, PosInterDist posInterDist, double epsilon) {

		var settings = new KStar.Settings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(1);

		if (generate) {

			settings.setEpsilon(this.epsilon);
			//settings.setShowPfuncProgress(true);

			settings.addScoreConsoleWriter(info -> {

				Function<PartitionFunction.Result,String> formatPfunc = (pfuncResult) -> {
					if (pfuncResult.status == PartitionFunction.Status.Estimated) {
						var bounds = pfuncResult.values.calcBounds();
						return String.format("\"%-12s\", \"%-12s\"",
							Log.formatBigEngineering(bounds.lower),
							Log.formatBigEngineering(bounds.upper)
						);
					}
					return "null";
				};

				return String.format("assertSequence(result, %3d, \"%s\", %s, %s, %s);",
					info.sequenceNumber,
					info.sequence.toString(Sequence.Renderer.ResType),
					formatPfunc.apply(info.kstarScore.protein),
					formatPfunc.apply(info.kstarScore.ligand),
					formatPfunc.apply(info.kstarScore.complex)
				);
			});
		}

		KStar kstar = new KStar(protein, ligand, complex, settings.build());

		var parallelism = Parallelism.makeCpu(4);
		try (var tasks = parallelism.makeTaskExecutor()) {

			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				ConfEnergyCalculator ecalc = ConfEnergyCalculator.makeBest(confSpace, parallelism);

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

		// the two pfunc bounds should overlap for each state
		assertThat(seq.score.protein.values.calcBounds().intersects(proteinExp), is(true));
		assertThat(seq.score.ligand.values.calcBounds().intersects(ligandExp), is(true));
		assertThat(seq.score.complex.values.calcBounds().intersects(complexExp), is(true));
	}

	public static void main(String[] args) {

		var test = new TestKStar();
		test.generate = true;

		test.epsilon = 0.1;
		test.test2RL0_TighterBounds();
	}
}
