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
import edu.duke.cs.osprey.kstar.*;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.List;


public class TestBBKStar {

	// import the compiled ccsx (from PINT or GUI) and set epsilon

	@Test
	public void test2RL0() {

		ConfSpace complex = ConfSpace.fromBytes(FileTools.readResourceBytes("/Thanatin/two-mut-D-peptide/complex.ccsx"));
		ConfSpace protein = ConfSpace.fromBytes(FileTools.readResourceBytes("/Thanatin/two-mut-D-peptide/CALP.ccsx"));
		ConfSpace ligand = ConfSpace.fromBytes(FileTools.readResourceBytes("/Thanatin/two-mut-D-peptide/match4.ccsx"));

		final double epsilon = 0.05;
		run(ligand, protein, complex, epsilon);
	}

	private static int run(ConfSpace protein, ConfSpace ligand, ConfSpace complex, double epsilon) {

		// format Kstar score information
		KStarScoreWriter.Formatter testFormatter = info ->
			String.format("%3d %s   protein: %s   ligand: %s   complex: %s   K*: %s",
				info.sequenceNumber,
				info.sequence.toString(Sequence.Renderer.ResType),
				info.kstarScore.protein.toString(),
				info.kstarScore.ligand.toString(),
				info.kstarScore.complex.toString(),
				info.kstarScore.toString()
			);

		// customize Kstar and BBKstar settings
		KStarSettings kstarSettings = new KStarSettings.Builder()
			.setEpsilon(epsilon)
			.setStabilityThreshold(null)
				// change this value for large designs
			.setMaxSimultaneousMutations(4)
				// set false to disable lots of printouts
			.setShowPfuncProgress(true)
			.addScoreConsoleWriter(testFormatter)
			.build();
		BBKStar.Settings bbkstarSettings = new BBKStar.Settings.Builder()
				// # of best seqs before BBK* stops
			.setNumBestSequences(4)
				// make this even multiple of available threads
			.setNumConfsPerBatch(8)
			.build();
		BBKStar bbkstar = new BBKStar(protein, ligand, complex, kstarSettings, bbkstarSettings);

		try (ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor()) {
			tasks.start(4);

			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				// turn off default confDB for tests - I have no idea what this does
				// info.confDBFile = null;

				// pass ConfSpace info?
				ConfSpace confSpace = (ConfSpace)info.confSpace;

				// determine how residue interactions distributed among fragments
				PosInterDist posInterDist = PosInterDist.DesmetEtAl1992;
				boolean includeStaticStatic = true;
				ConfEnergyCalculator ecalc = new CPUConfEnergyCalculator(confSpace);

				// calculate minimized reference energies
				SimpleReferenceEnergies eref = new ErefCalculator.Builder(ecalc)
					.setMinimize(true)
					.build()
					.calc();

				// create energy calculator amenable to ccsx file format
				info.confEcalcMinimized = new ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)
					.setPosInterDist(posInterDist)
					.setReferenceEnergies(eref)
					.setMinimize(true)
					.setIncludeStaticStatic(includeStaticStatic)
					.build();

				// build energy matrix + A* search tree
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

				// use gradient descent to min pfuncs
				info.pfuncFactory = rcs -> new GradientDescentPfunc(
					info.confEcalcMinimized,
					info.confSearchFactoryMinimized.make(rcs),
					info.confSearchFactoryMinimized.make(rcs),
					rcs.getNumConformations()
				).setPreciseBcalc(true);
			}


			// run K*
			List<ScoredSequence> sequences = bbkstar.run(tasks);
			int lsize = sequences.size();
			System.out.println("Length of sequences: " + lsize);

			// make Sequence Analyzer + save PDB ensembles
			SequenceAnalyzer analyzer = new SequenceAnalyzer(bbkstar);
			for (ScoredSequence sequence : sequences) {
				System.out.println("result:");
				System.out.println("\tsequence: " + sequence.sequence());
				System.out.println("\tscore: " + sequence.score());

				// set # conformations printed in ensemble + analyze
				int numEConfs = 10;
				SequenceAnalyzer.Analysis analysis = analyzer.analyze(sequence.sequence(), numEConfs);
				// System.out.println(analysis);

				// formats seqstr for file outputs (only changes filename)
				String seqstr = sequence.sequence().toString(Sequence.Renderer.ResTypeMutations)
						.replace(' ', '-');

				// set where file is saved + name
				File ensembleFile = new File(String.format("/home/henry-childs/IdeaProjects/OSPREY3/src/test/resources/Thanatin/two-mut-D-peptide/seq.%s.pdb", seqstr));

				// write the PDB (can also set filepath here)
				analysis.writePdb(ensembleFile.getAbsolutePath(), String.format("Top %d conformations for sequence %s",
						numEConfs, sequence.sequence()));
			}

			return 1;

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
