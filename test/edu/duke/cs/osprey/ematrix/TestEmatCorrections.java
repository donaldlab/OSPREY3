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

package edu.duke.cs.osprey.ematrix;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.Test;

import java.util.*;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static edu.duke.cs.osprey.tools.Log.log;


public class TestEmatCorrections {

	// conf spaces with < 3 positions won't get corrected, but should still be accurate
	@Test
	public void testAccurateTriples_1CC8_2N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional,
			3
		);
	}
	@Test
	public void testAccurateQuads_1CC8_2N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional,
			4
		);
	}

	// AllOnPairs conf spaces with 2 positions won't get corrected either, but the bounds should be perfect anyway! =F
	@Test
	public void testPerfectTriples_1CC8_2N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}
	@Test
	public void testPerfectTriples_1CC8_2F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}
	@Test
	public void testPerfectQuads_1CC8_2N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}
	@Test
	public void testPerfectQuads_1CC8_2F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}

	// conf spaces with >= 3 positions should get corrected by triples, and the corrections should be accurate
	@Test
	public void testAccurateTriples_1CC8_3N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional,
			3
		);
	}
	@Test
	public void testAccurateTriples_1CC8_3F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.Traditional,
			3
		);
	}

	// conf spaces with >= 4 positions should get corrected by quads, and the corrections should be accurate
	@Test
	public void testAccurateQuads_1CC8_4N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional,
			4
		);
	}
	@Test
	public void testAccurateQuads_1CC8_4F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional,
			4
		);
	}

	// triples corrections should be perfect for AllOnPairs with 3 positions! =D
	@Test
	public void testPerfectTriples_1CC8_3N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}
	@Test
	public void testPerfectTriples_1CC8_3F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47" // val
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}

	// quads corrections should be perfect for AllOnPairs with 4 positions! =D
	@Test
	public void testPerfectQuads_1CC8_4N_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}
	@Test
	public void testPerfectQuads_1CC8_4F_AllOnPairs() {
		assertPerfectCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}

	// triples should help for n > 3 residues
	@Test
	public void testAccurateTriples_1CC8_4N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional,
			3
		);
	}
	@Test
	public void testAccurateTriples_1CC8_4F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.Traditional,
			3
		);
	}
	@Test
	public void testAccurateTriples_1CC8_4N_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}
	@Test
	public void testAccurateTriples_1CC8_4F_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60" // ile
			),
			EnergyPartition.AllOnPairs,
			3
		);
	}

	// quads should help for n > 4 residues
	@Test
	public void testAccurateQuads_1CC8_5N_Traditional() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60", // ile
				"A56" // ile
			),
			EnergyPartition.Traditional,
			4
		);
	}
	@Test
	public void testAccurateQuads_1CC8_5F_Traditional() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60", // ile
				"A56" // ile
			),
			EnergyPartition.Traditional,
			4
		);
	}
	@Test
	public void testAccurateQuads_1CC8_5N_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpaceNoFixed(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60", // ile
				"A56" // ile
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}
	@Test
	public void testAccurateQuads_1CC8_5F_AllOnPairs() {
		assertAccurateCorrections(
			makeConfSpace(
				"/1CC8.ss.pdb",
				"A26", // leu
				"A45", // val
				"A47", // val
				"A60", // ile
				"A56" // ile
			),
			EnergyPartition.AllOnPairs,
			4
		);
	}


	private static SimpleConfSpace makeConfSpace(String pdbPath, String ... resNums) {

		// share the template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();

		// build strands, one for each residue
		// so we don't end up with any fixed residues
		Molecule pdb = PDBIO.readResource(pdbPath);
		Strand strand = new Strand.Builder(pdb)
				.setTemplateLibrary(templateLib)
				.build();
		for (String resNum : resNums) {
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		return new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
	}

	private static SimpleConfSpace makeConfSpaceNoFixed(String pdbPath, String ... resNums) {

		// share the template library
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder().build();

		// build strands, one for each residue
		// so we don't end up with any fixed residues
		Molecule pdb = PDBIO.readResource(pdbPath);
		List<Strand> strands = new ArrayList<>();
		for (String resNum : resNums) {
			Strand strand = new Strand.Builder(pdb)
				.setTemplateLibrary(templateLib)
				.setResidues(resNum, resNum)
				.build();
			strands.add(strand);
			strand.flexibility.get(resNum)
				.setLibraryRotamers(Strand.WildType)
				.addWildTypeRotamers()
				.setContinuous();
		}

		return new SimpleConfSpace.Builder()
			.addStrands(strands)
			.build();
	}

	private static void assertPerfectCorrections(SimpleConfSpace confSpace, EnergyPartition epart, int order) {
		forEachTopConf(confSpace, epart, order, (index, conf, energy, lowerBound, correctedBound) -> {
			assertThat(energy - correctedBound, isAbsolutely(0, 1e-4));
		});
	}

	private static void assertAccurateCorrections(SimpleConfSpace confSpace, EnergyPartition epart, int order) {
		forEachTopConf(confSpace, epart, order, (index, conf, energy, lowerBound, correctedBound) -> {
			assertThat(energy - correctedBound, greaterThanOrEqualTo(-1e-4));
		});
	}

	interface ConfListener {
		void onConf(int index, int[] conf, double energy, double lowerBound, double correctedBound);
	}

	private static void forEachTopConf(SimpleConfSpace confSpace, EnergyPartition epart, int order, ConfListener listener) {

		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.setParallelism(Parallelism.makeCpu(4))
			.build()) {

			ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
				.setEnergyPartition(epart)
				.build();

			// calc an emat
			EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setTripleCorrectionThreshold(order >= 3 ? Double.POSITIVE_INFINITY : null)
				.setQuadCorrectionThreshold(order >= 4 ? Double.POSITIVE_INFINITY : null)
				.build()
				.calcEnergyMatrix();

			// get the top 10 confs by energy lower bound
			ConfAStarTree astar = new ConfAStarTree.Builder(emat, confSpace)
				.setTraditional()
				.build();
			for (int i=0; i<10; i++) {

				// get the conf and the lower bound
				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					break;
				}

				// calculate the energy
				ConfSearch.EnergiedConf econf = confEcalc.calcEnergy(conf);

				// calculate the corrected bound
				double[] corrected = { econf.getScore() };
				emat.forEachHigherOrderTupleIn(econf.getAssignments(), (tuple, tupleEnergy) -> {
					corrected[0] += tupleEnergy;
				});

				// TEMP
				log("conf %2d   e=%9.3f   lb=%9.3f (%6.3f)   clb=%9.3f (%6.3f)",
					i,
					econf.getEnergy(),
					econf.getScore(), econf.getEnergy() - econf.getScore(),
					corrected[0], econf.getEnergy() - corrected[0]
				);

				listener.onConf(
					i,
					econf.getAssignments(),
					econf.getEnergy(),
					econf.getScore(),
					corrected[0]
				);
			}
		}
	}
}
