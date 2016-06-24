package edu.duke.cs.osprey;

import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;

public class TestBase {
	
	private static final double Epsilon = 1e-14;

	protected double getRelativeError(double expected, double observed) {
		return Math.abs(expected - observed)/Math.abs(observed);
	}
	
	protected Matcher<Double> isRelatively(final double expected) {
		return new BaseMatcher<Double>() {

			@Override
			public boolean matches(Object obj) {
				double observed = ((Double)obj).doubleValue();
				return getRelativeError(expected, observed) <= Epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				double observed = ((Double)obj).doubleValue();
				double relErr = getRelativeError(expected, observed);
				desc.appendText("value ").appendValue(observed)
					.appendText(" has relative err ").appendValue(relErr)
					.appendText(" that's greater than epsilon ").appendValue(Epsilon);
			}
		};
	}
	
	protected static ForcefieldParams makeDefaultFFParams() {
		// values from default config file
		String forceField = "AMBER";
		boolean distDeptDielect = true;
		double dielectConst = 6;
		double vdwMult = 0.95;
		boolean doSolv = true;
		double solvScale = 0.5;
		boolean useHForElectrostatics = true;
		boolean useHForVdw = true;
		return new ForcefieldParams(
			forceField, distDeptDielect, dielectConst, vdwMult,
			doSolv, solvScale, useHForElectrostatics, useHForVdw
		);
	}
	
	protected static void initDefaultEnvironment() {
		
		// values from default config file
		EnvironmentVars.setDataDir("dataFiles");
		
		// make energy function
		double shellDistCutoff = Double.POSITIVE_INFINITY;
		boolean usePoissonBoltzmann = false;
		EnvironmentVars.curEFcnGenerator = new EnergyFunctionGenerator(makeDefaultFFParams(), shellDistCutoff, usePoissonBoltzmann);
		
		// make residue templates
		EnvironmentVars.resTemplates = new GenericResidueTemplateLibrary(
			new String[] { "all_amino94.in", "all_aminont94.in", "all_aminoct94.in", "all_nuc94_and_gr.in" },
			makeDefaultFFParams()
		);
		EnvironmentVars.resTemplates.loadTemplateCoords("all_amino_coords.in");
		EnvironmentVars.resTemplates.makeDAminoAcidTemplates();
		
		// make rotamers
		boolean useBackboneDependentRotamers = false;
		EnvironmentVars.resTemplates.loadRotamerLibrary("LovellRotamer.dat", useBackboneDependentRotamers);
	}
}
