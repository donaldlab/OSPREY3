package edu.duke.cs.osprey.kstar.compiled;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
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
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.Timeout;

import java.util.List;
import java.util.concurrent.TimeUnit;


public class TestKStar {

	@Rule
	public Timeout globalTimeout = new Timeout(5, TimeUnit.MINUTES);

	// NOTE: these tests don't test for correctness, only that the code finishes and doesn't crash

	@Test
	public void test2RL0() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"));

		final double epsilon = 0.99;
		run(chainG, chainA, complex, epsilon);
	}

	@Test
	public void testF98Y() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.complex.ccsx"));
		ConfSpace dhfr = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.DHFR.ccsx"));
		ConfSpace nadph06w = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.NADPH.06W.ccsx"));

		final double epsilon = 0.05;
		run(dhfr, nadph06w, complex, epsilon);
	}

	private static List<KStar.ScoredSequence> run(ConfSpace protein, ConfSpace ligand, ConfSpace complex, double epsilon) {

		KStarScoreWriter.Formatter testFormatter = info ->
			String.format("%3d %s   protein: %s   ligand: %s   complex: %s   K*: %s",
				info.sequenceNumber,
				info.sequence.toString(Sequence.Renderer.ResType),
				info.kstarScore.protein.toString(),
				info.kstarScore.ligand.toString(),
				info.kstarScore.complex.toString(),
				info.kstarScore.toString()
			);

		KStar.Settings settings = new KStar.Settings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(1)
			//.setShowPfuncProgress(true)
			.addScoreConsoleWriter(testFormatter)
			.build();
		KStar kstar = new KStar(protein, ligand, complex, settings);

		try (ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(4);

			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				// turn off the default confdb for tests
				info.confDBFile = null;

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
				boolean minimize = true;
				boolean includeStaticStatic = true;
				ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

				SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
					.setMinimize(minimize)
					.build()
					.calc();

				EnergyMatrix emat = new EmatCalculator.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(minimize)
					.setIncludeStaticStatic(includeStaticStatic)
					.build()
					.calc();

				info.confEcalc = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(minimize)
					.setIncludeStaticStatic(includeStaticStatic)
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
			return kstar.run(tasks);

		} finally {

			// cleanup
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {
				if (info.confEcalc != null) {
					((ConfEnergyCalculatorAdapter)info.confEcalc).confEcalc.close();
				}
			}
		}
	}
}
