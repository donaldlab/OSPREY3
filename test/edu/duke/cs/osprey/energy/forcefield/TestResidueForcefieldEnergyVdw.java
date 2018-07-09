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

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Residues;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.ArrayList;
import java.util.function.Function;

public class TestResidueForcefieldEnergyVdw {

	@BeforeClass
	public static void before() {
		TestForcefieldEnergy.before();
	}

	public static void main(String[] args) {

		before();

		TestForcefieldEnergy.TestResidues r = new TestForcefieldEnergy.TestResidues();
		Residues residues = new Residues(
			//r.gly15
			//r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		);

		new EnergyCalculator.Builder(residues, new ForcefieldParams())
			.setType(EnergyCalculator.Type.Cpu)
			.setParallelism(Parallelism.makeCpu(1))
			.use((ecalc) -> {

				// make complete interactions
				ResidueInteractions inters = new ResidueInteractions();
				inters.addComplete(residues);

				// get the energy function
				ParametricMolecule pmol = new ParametricMolecule(r.mol, new ArrayList<>(), null);
				ResidueForcefieldEnergy efunc = (ResidueForcefieldEnergy)ecalc.makeEnergyFunction(pmol, inters);

				System.out.println(String.format("vdw: %20.12f", efunc.getVanDerWaalsEnergy()));
			});
	}

	private void check(Function<TestForcefieldEnergy.TestResidues,Residues> residuesFactory, double expected) {

		Residues residues = residuesFactory.apply(new TestForcefieldEnergy.TestResidues());

		ResPairCache resPairCache = new ResPairCache(
			new ForcefieldParams(),
			new AtomConnectivity.Builder()
				.addTemplates(residues)
				.build()
		);
		ResidueInteractions inters = new ResidueInteractions();
		inters.addComplete(residues);

		ResidueForcefieldEnergy.Vdw vdw = new ResidueForcefieldEnergy.Vdw(resPairCache, inters, residues);
		double energy = vdw.getEnergy();

		assertThat(energy, isAbsolutely(expected, 1e-12));
	}

	@Test
	public void singleGly() {
		check((r) -> new Residues(
			r.gly15
		), -0.122472744448);
	}

	@Test
	public void glyPair() {
		check((r) -> new Residues(
			r.gly06, r.gly15
		), -0.260879089375);
	}

	@Test
	public void glySerPair() {
		check((r) -> new Residues(
			r.gly15, r.ser17
		), 0.810077899104);
	}

	@Test
	public void trpPair() {
		check((r) -> new Residues(
			r.trp18, r.trp25
		), -0.002972265815);
	}

	@Test
	public void the4Residues() {
		check((r) -> new Residues(
			r.gly15, r.ser17, r.trp18, r.trp25
		), -1.931665689920);
	}

	@Test
	public void the6Residues() {
		check((r) -> new Residues(
			r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24
		), -7.805339930453);
	}

	@Test
	public void the10Residues() {
		check((r) -> new Residues(
			r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34
		), -9.988696541744);
	}

	@Test
	public void the14Residues() {
		check((r) -> new Residues(
			r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48
		), -13.153249044828);
	}

	@Test
	public void the24Residues() {
		check((r) -> new Residues(
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		), -12.562797940805);
	}
}
