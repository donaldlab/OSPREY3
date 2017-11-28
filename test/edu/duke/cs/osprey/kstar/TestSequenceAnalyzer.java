package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.Test;

public class TestSequenceAnalyzer {

	private static final double EnergyEpsilon = 1e-6;

	@Test
	public void test2RL0() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		Parallelism parallelism = Parallelism.makeCpu(4);
		//Parallelism parallelism = Parallelism.make(4, 1, 8);

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
				KStar.ConfSearchFactory confSearchFactory = (emat, pmat) -> {
					return new ConfAStarTree.Builder(emat, pmat)
						.setTraditional()
						.build();
				};

				KStar.Settings settings = new KStar.Settings.Builder().build();

				SequenceAnalyzer analyzer = new SequenceAnalyzer(
					confSpaces.protein,
					confSpaces.ligand,
					confSpaces.complex,
					ecalc,
					confEcalcFactory,
					confSearchFactory,
					settings
				);

				SequenceAnalyzer.Analysis analysis = analyzer.analyze(
					confSpaces.complex.makeWildTypeSequence(),
					1
				);

				assertThat(analysis.info.type, is(KStar.ConfSpaceType.Complex));
				assertThat(analysis.sequence.toString(Sequence.Renderer.ResTypeMutations), is("phe asp glu thr phe lys ile thr"));
				assertThat(analysis.ensemble.analyses.size(), is(4));
				assertThat(analysis.ensemble.analyses.get(0).epmol.energy, isAbsolutely(-68.416267, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(1).epmol.energy, isAbsolutely(-68.244583, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(2).epmol.energy, isAbsolutely(-68.114214, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(3).epmol.energy, isAbsolutely(-67.945550, EnergyEpsilon));

				analysis = analyzer.analyze(
					confSpaces.complex.makeWildTypeSequence().set("G649", "ILE"),
					1
				);

				assertThat(analysis.info.type, is(KStar.ConfSpaceType.Complex));
				assertThat(analysis.sequence.toString(Sequence.Renderer.ResTypeMutations), is("ILE asp glu thr phe lys ile thr"));
				assertThat(analysis.ensemble.analyses.size(), is(12));
				assertThat(analysis.ensemble.analyses.get(0).epmol.energy, isAbsolutely(-61.938135, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(1).epmol.energy, isAbsolutely(-61.905580, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(2).epmol.energy, isAbsolutely(-61.890854, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(3).epmol.energy, isAbsolutely(-61.639204, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(4).epmol.energy, isAbsolutely(-61.637555, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(5).epmol.energy, isAbsolutely(-61.608625, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(6).epmol.energy, isAbsolutely(-61.593341, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(7).epmol.energy, isAbsolutely(-61.591667, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(8).epmol.energy, isAbsolutely(-61.586393, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(9).epmol.energy, isAbsolutely(-61.335184, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(10).epmol.energy, isAbsolutely(-61.311892, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(11).epmol.energy, isAbsolutely(-61.247090, EnergyEpsilon));

				analysis = analyzer.analyze(
					confSpaces.protein.makeWildTypeSequence().set("G649", "ILE"),
					1
				);

				assertThat(analysis.info.type, is(KStar.ConfSpaceType.Protein));
				assertThat(analysis.sequence.toString(Sequence.Renderer.ResTypeMutations), is("ILE asp glu thr"));
				assertThat(analysis.ensemble.analyses.size(), is(6));
				assertThat(analysis.ensemble.analyses.get(0).epmol.energy, isAbsolutely(-3.000911, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(1).epmol.energy, isAbsolutely(-2.903385, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(2).epmol.energy, isAbsolutely(-2.508572, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(3).epmol.energy, isAbsolutely(-2.418352, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(4).epmol.energy, isAbsolutely(-2.401930, EnergyEpsilon));
				assertThat(analysis.ensemble.analyses.get(5).epmol.energy, isAbsolutely(-2.309892, EnergyEpsilon));
			});
	}
}
