/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.kstar;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;
import static edu.duke.cs.osprey.TestBase.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
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
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
			.setParallelism(parallelism)
			.build()) {

			KStar.Settings settings = new KStar.Settings.Builder().build();
			KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				info.confEcalc = new ConfEnergyCalculator.Builder(info.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(info.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					).build();

				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
					.build()
					.calcEnergyMatrix();

				info.confSearchFactory = (rcs) ->
					new ConfAStarTree.Builder(emat, rcs)
						.setTraditional()
						.build();
			}

			SequenceAnalyzer analyzer = new SequenceAnalyzer(kstar);

			SequenceAnalyzer.Analysis analysis = analyzer.analyze(
				confSpaces.complex.makeWildTypeSequence(),
				1
			);

			assertThat(analysis.info.id, is("complex"));
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

			assertThat(analysis.info.id, is("complex"));
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

			assertThat(analysis.info.id, is("protein"));
			assertThat(analysis.sequence.toString(Sequence.Renderer.ResTypeMutations), is("ILE asp glu thr"));
			assertThat(analysis.ensemble.analyses.size(), is(6));
			assertThat(analysis.ensemble.analyses.get(0).epmol.energy, isAbsolutely(-3.000911, EnergyEpsilon));
			assertThat(analysis.ensemble.analyses.get(1).epmol.energy, isAbsolutely(-2.903385, EnergyEpsilon));
			assertThat(analysis.ensemble.analyses.get(2).epmol.energy, isAbsolutely(-2.508572, EnergyEpsilon));
			assertThat(analysis.ensemble.analyses.get(3).epmol.energy, isAbsolutely(-2.418352, EnergyEpsilon));
			assertThat(analysis.ensemble.analyses.get(4).epmol.energy, isAbsolutely(-2.401930, EnergyEpsilon));
			assertThat(analysis.ensemble.analyses.get(5).epmol.energy, isAbsolutely(-2.309892, EnergyEpsilon));
		}
	}
}
