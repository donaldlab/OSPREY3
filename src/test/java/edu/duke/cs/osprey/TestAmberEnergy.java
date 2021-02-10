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

package edu.duke.cs.osprey;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.structure.PDBIO;

public class TestAmberEnergy extends TestBase {
	
	// NOTE: for tests, try to rely on as little Dynamic global configuration as possible
	// Static config (stuff that never changes between runs) is ok
	// we want to isolate the test case as much as possible from the rest of the codebase
	// ideally, the test should completely configure itself and rely on no external settings
	// the test should be capable of running correctly independently without needing config from the user
	// ie, don't rely on values from the config file parser!
	
	// also, don't write to logs
	// use assert statements to check for correct output
	
	@Test
	public void test1CC8WithSolv() {
		test("examples/1CC8/1CC8.ss.pdb", SolvationForcefield.EEF1, -986.6491862981836);
	}

	@Test
	public void test1CC8NoSolv() {
		test("examples/1CC8/1CC8.ss.pdb", null, -639.7025085949941);
	}

	private void test(String pdbPath, SolvationForcefield solvff, double energy) {
		
		//compute the full energy for 1CC8 using the default AMBER forcefield, and compare it to OSPREY 2 values
		Strand strand = new Strand.Builder(PDBIO.readFile(pdbPath)).build();
		
		ForcefieldParams ffparams = new ForcefieldParams();
		ffparams.solvationForcefield = solvff;
		EnergyFunction efunc = new EnergyFunctionGenerator(ffparams).fullMolecEnergy(strand.mol);
		
		assertThat(efunc.getEnergy(), isRelatively(energy));
	}
}
