package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;

public class TestAmberEnergy extends TestBase {
	
	// NOTE: for tests, try to rely on as little Dynamic global configuration as possible
	// Static config (stuff that never changes between runs) is ok
	// we want to isolate the test case as much as possible from the rest of the codebase
	// ideally, the test should completely configure itself and rely on no external settings
	// the test should be capable of running correctly independently without needing config from the user
	// ie, don't rely on values from the config file parser!
	
	// also, don't write to logs
	// use assert statements to check for correct output
	
	@BeforeClass
	public static void before() {
            
		// load static config for amber forcefield
		EnvironmentVars.setDataDir("dataFiles");
		EnvironmentVars.resTemplates = new GenericResidueTemplateLibrary(
			new String[] { "all_amino94.in", "all_aminont94.in", "all_aminoct94.in", "all_nuc94_and_gr.in" },
			makeFFParams(false)
		);
	}
	
	private static ForcefieldParams makeFFParams(boolean doSolv) {
		// NOTE: values from default config file
		String forceField = "AMBER";
		boolean distDeptDielect = true;
		double dielectConst = 6;
		double vdwMult = 0.95;
		double solvScale = 0.5;
		boolean useHForElectrostatics = true;
		boolean useHForVdw = true;
		return new ForcefieldParams(
			forceField, distDeptDielect, dielectConst, vdwMult,
			doSolv, solvScale, useHForElectrostatics, useHForVdw
		);
	}
	
        @Test
	public void test1CC8WithSolv()
	throws Exception {
		
		// setup the test
		Molecule m = PDBFileReader.readPDBFile("test/1CC8/1CC8.ss.pdb");
		assertThat(m, is(not(nullValue())));
		EnergyFunction efunc = makeEfunc(m, true);
		
		// run the test
		double energy = efunc.getEnergy();
		
		// check the result
		assertThat(energy, isRelatively(-986.6491862981836));
	}

	@Test
	public void test1CC8NoSolv()
	throws Exception {
		
		// setup the test
		Molecule m = PDBFileReader.readPDBFile("test/1CC8/1CC8.ss.pdb");
		assertThat(m, is(not(nullValue())));
		EnergyFunction efunc = makeEfunc(m, false);
		
		// run the test
		double energy = efunc.getEnergy();
		
		// check the result
		assertThat(energy, isRelatively(-639.7025085949941));
	}
	
	private EnergyFunction makeEfunc(Molecule mol, boolean doSolv) {
		double shellDistCutoff = Double.POSITIVE_INFINITY;
		boolean usePoissonBoltzmann = false;
		EnergyFunctionGenerator efuncgen = new EnergyFunctionGenerator(
			makeFFParams(doSolv),
			shellDistCutoff,
			usePoissonBoltzmann
		);
		return efuncgen.fullMolecEnergy(mol);
	}
}
