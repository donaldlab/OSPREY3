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
