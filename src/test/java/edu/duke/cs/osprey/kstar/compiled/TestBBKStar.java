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
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
//import org.junit.Rule;
import org.junit.jupiter.api.Test;
//import org.junit.rules.Timeout;

import java.util.List;
import java.util.concurrent.TimeUnit;

public class TestBBKStar {

	// NOTE: these tests don't test for correctness, only that the code finishes and doesn't crash

	@Test
	public void test2RL0() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccsx"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccsx"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccsx"));

		final double epsilon = 0.99;
		run(chainG, chainA, complex, epsilon);
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

		KStar.Settings kstarSettings = new KStar.Settings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
			.setMaxSimultaneousMutations(1)
			//.setShowPfuncProgress(true)
			.addScoreConsoleWriter(testFormatter)
			.build();
		BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
			.setNumBestSequences(5)
			.setNumConfsPerBatch(8)
			.build();
		BBKStar bbkstar = new BBKStar(protein, ligand, complex, kstarSettings, bbkstarSettings);

		try (ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(4);

			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				// turn off default confDB for tests
				info.confDBFile = null;

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
				boolean includeStaticStatic = true;
				ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

				SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
					.setMinimize(true)
					.build()
					.calc();

				info.confEcalcMinimized = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(true)
					.setIncludeStaticStatic(includeStaticStatic)
					.build();

				EnergyMatrix ematMinimized = new EmatCalculator.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(true)
					.setIncludeStaticStatic(includeStaticStatic)
					.build()
					.calc();
				info.confSearchFactoryMinimized = (rcs) ->
					new ConfAStarTree.Builder(ematMinimized, rcs)
						.setTraditional()
						.build();

				// BBK* needs rigid energies too
				EnergyMatrix ematRigid = new EmatCalculator.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(false)
					.setIncludeStaticStatic(includeStaticStatic)
					.build()
					.calc();
				info.confSearchFactoryRigid = (rcs) ->
					new ConfAStarTree.Builder(ematRigid, rcs)
						.setTraditional()
						.build();

				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalcMinimized,
					info.confSearchFactoryMinimized.make(rcs),
					info.confSearchFactoryMinimized.make(rcs),
					rcs.getNumConformations()
				).setPreciseBcalc(true);
			}

			// run K*
			return bbkstar.run(tasks);

		} finally {

			// cleanup
			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {
				if (info.confEcalcMinimized != null) {
					((ConfEnergyCalculatorAdapter)info.confEcalcMinimized).confEcalc.close();
				}
			}
		}
	}
}
