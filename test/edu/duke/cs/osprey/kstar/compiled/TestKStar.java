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
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.util.List;


public class TestKStar {

	@Test
	public void test2RL0() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.complex.ccs.xz"));
		ConfSpace chainA = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.A.ccs.xz"));
		ConfSpace chainG = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/2RL0.G.ccs.xz"));

		final double epsilon = 0.99;
		run(chainG, chainA, complex, epsilon);
	}

	@Test
	public void testF98Y() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.complex.ccs.xz"));
		ConfSpace dhfr = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.DHFR.ccs.xz"));
		ConfSpace nadph06w = ConfSpace.fromBytes(FileTools.readFileBytes("examples/python.ccs/F98Y/4tu5.NADPH.06W.ccs.xz"));

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

				ConfSpace confSpace = (ConfSpace)info.confSpace;

				PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
				boolean minimize = true;
				boolean includeStaticStatic = true;
				ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace, tasks);

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

				info.confEcalc = new ConfEnergyCalculatorAdapter.Builder(ecalc)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(minimize)
					.setIncludeStaticStatic(includeStaticStatic)
					.build();

				info.confSearchFactory = (rcs) ->
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();
			}

			// run K*
			return kstar.run();

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
