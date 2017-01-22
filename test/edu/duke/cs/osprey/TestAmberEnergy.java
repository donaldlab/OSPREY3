package edu.duke.cs.osprey;

import static org.junit.Assert.*;

import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
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
	public void test1CC8WithSolv()
	throws Exception {
		
		// setup the test
		Strand strand = Strand.builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		EnergyFunction efunc = makeEfunc(strand, true);
		
		// run the test
		double energy = efunc.getEnergy();
		
		// check the result
		assertThat(energy, isRelatively(-986.6491862981836));
	}

	@Test
	public void test1CC8NoSolv()
	throws Exception {
		
		// setup the test
		Strand strand = Strand.builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build();
		EnergyFunction efunc = makeEfunc(strand, false);
		
		// run the test
		double energy = efunc.getEnergy();
		
		// check the result
		assertThat(energy, isRelatively(-639.7025085949941));
	}
	
	private EnergyFunction makeEfunc(Strand strand, boolean doSolv) {
		ForcefieldParams ffparams = new ForcefieldParams(strand.templateLib.ffParams);
		ffparams.doSolvationE = doSolv;
		// TODO: move these params into ForcefieldParams
		double shellDistCutoff = Double.POSITIVE_INFINITY;
		boolean usePoissonBoltzmann = false;
		EnergyFunctionGenerator efuncgen = new EnergyFunctionGenerator(
			ffparams,
			shellDistCutoff,
			usePoissonBoltzmann
		);
		return efuncgen.fullMolecEnergy(strand.mol);
	}
}
