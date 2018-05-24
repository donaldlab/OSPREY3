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

package edu.duke.cs.osprey.energy.forcefield;

import static edu.duke.cs.osprey.TestBase.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Residues;
import org.junit.BeforeClass;
import org.junit.Test;

public class TestResidueForcefieldSubEnergies {

	@BeforeClass
	public static void before() {
		TestForcefieldEnergy.before();
	}

	private static void checkEnergies(Residues residues) {

		ResidueInteractions inters = TestForcefieldEnergy.IntersType.AllPairs.makeInters(residues);

		for (TestForcefieldEnergy.FFType fftype : TestForcefieldEnergy.FFType.values()) {

			// build the energy function
			ForcefieldParams ffparams = fftype.makeFFParams();
			AtomConnectivity connectivity = new AtomConnectivity.Builder()
				.addTemplates(residues)
				.setParallelism(Parallelism.makeCpu(4))
				.build();
			ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
			ResidueForcefieldEnergy efunc = new ResidueForcefieldEnergy(resPairCache, inters, residues);

			// add up the component energies
			double totalSum = 0.0;
			for (ResidueInteractions.Pair resPair : inters) {

				ResidueForcefieldEnergy pairEfunc = new ResidueForcefieldEnergy(resPairCache, new ResidueInteractions(resPair), residues);

				double pairSum = 0.0;
				pairSum += pairEfunc.getElectrostaticsEnergy();
				pairSum += pairEfunc.getVanDerWaalsEnergy();
				pairSum += pairEfunc.getSolvationEnergy();
				pairSum += pairEfunc.getOffsetsEnergy();

				double pairEnergy = pairEfunc.getEnergy();

				assertThat("res pair: " + resPair.resNum1 + ":" + resPair.resNum2 + ", forcefield type: " + fftype,
					pairSum, isAbsolutely(pairEnergy, 1e-12)
				);

				totalSum += pairSum;
			}

			double totalEnergy = efunc.getEnergy();
			assertThat("forcefield type: " + fftype,
				totalSum, isAbsolutely(totalEnergy, 1e-12)
			);
		}
	}

	@Test
	public void singleGly() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15));
	}

	@Test
	public void glyPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly06, r.gly15));
	}

	@Test
	public void glySerPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15, r.ser17));
	}

	@Test
	public void trpPair() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.trp18, r.trp25));
	}

	@Test
	public void the4Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15, r.ser17, r.trp18, r.trp25));
	}

	@Test
	public void the6Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24));
	}

	@Test
	public void the10Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34));
	}

	@Test
	public void the14Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48));
	}

	@Test
	public void the24Residues() {
		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		checkEnergies(new Residues(
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		));
	}
}
