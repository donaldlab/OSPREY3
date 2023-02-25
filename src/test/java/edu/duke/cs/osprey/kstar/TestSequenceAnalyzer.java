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
import static org.hamcrest.MatcherAssert.*;
import static edu.duke.cs.osprey.TestBase.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.util.function.Function;

public class TestSequenceAnalyzer {

	private static final double EnergyEpsilon = 1e-6;

	@Test
	public void test2RL0() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();

		Parallelism parallelism = Parallelism.makeCpu(4);
		//Parallelism parallelism = Parallelism.make(4, 1, 8);

		final String confdbPattern = "kstar.%s.conf.db";
		try (TempFile proteinDBFile = new TempFile("kstar.protein.conf.db")) {
			try (TempFile ligandDBFile = new TempFile("kstar.ligand.conf.db")) {
				try (TempFile complexDBFile = new TempFile("kstar.complex.conf.db")) {

					// how should we compute energies of molecules?
					try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
						.setParallelism(parallelism)
						.build()) {

						KStar.Settings settings = new KStar.Settings.Builder()
							.setShowPfuncProgress(true)
							.build();
						KStar kstar = new KStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, settings);
						for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

							SimpleConfSpace confSpace = (SimpleConfSpace)info.confSpace;

							info.confDBFile = new File(String.format(confdbPattern, info.type.name().toLowerCase()));

							info.confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
								.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
									.build()
									.calcReferenceEnergies()
								).build();

							EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(info.confEcalc)
								.build()
								.calcEnergyMatrix();

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

						// score just a few sequences
						Sequence[] seqs = {
							confSpaces.complex.makeWildTypeSequence(),
							confSpaces.complex.makeWildTypeSequence().set("G649", "ILE"),
							confSpaces.protein.makeWildTypeSequence().set("G649", "ILE")
						};

						kstar.score(seqs[0], ecalc.tasks);
						kstar.score(seqs[1], ecalc.tasks);
						// don't score seq 2, it's just for the protein conf space

						assertThat(proteinDBFile.exists(), is(true));
						assertThat(ligandDBFile.exists(), is(true));
						assertThat(complexDBFile.exists(), is(true));

						SequenceAnalyzer analyzer = new SequenceAnalyzer(kstar);
						SequenceAnalyzer.Analysis analysis;

						// analyze seq 0
						analysis = analyzer.analyze(seqs[0], 4);
						assertThat(analysis.info.id, is("complex"));
						assertThat(analysis.sequence.toString(Sequence.Renderer.ResTypeMutations), is("phe asp glu thr phe lys ile thr"));
						assertThat(analysis.ensemble.analyses.size(), is(4));
						assertThat(analysis.ensemble.analyses.get(0).epmol.energy, isAbsolutely(-68.416267, EnergyEpsilon));
						assertThat(analysis.ensemble.analyses.get(1).epmol.energy, isAbsolutely(-68.244583, EnergyEpsilon));
						assertThat(analysis.ensemble.analyses.get(2).epmol.energy, isAbsolutely(-68.114214, EnergyEpsilon));
						assertThat(analysis.ensemble.analyses.get(3).epmol.energy, isAbsolutely(-67.945550, EnergyEpsilon));

						analysis = analyzer.analyze(seqs[1], 12);
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

						analysis = analyzer.analyze(seqs[2], 6);
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
		}
	}
}
