package edu.duke.cs.osprey;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestBase {
	
	private static final double DefaultEpsilon = 1e-14;
	
	private static Map<EnergyMatrixConfig,EnergyMatrix> m_energyMatrixCache;
	
	static {
		m_energyMatrixCache = new HashMap<>();
	}
	
	public static class EnergyMatrixConfig {
		
		public String pdbPath;
		public int numFlexible;
		public boolean addWtRots;
		public boolean doMinimize;
		
		@Override
		public boolean equals(Object other) {
			if (other instanceof EnergyMatrixConfig) {
				return equals((EnergyMatrixConfig)other);
			}
			return false;
		}
		
		public boolean equals(EnergyMatrixConfig other) {
			return this.pdbPath.equals(other.pdbPath)
				&& this.numFlexible == other.numFlexible
				&& this.addWtRots == other.addWtRots
				&& this.doMinimize == other.doMinimize;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				this.pdbPath.hashCode(),
				Integer.valueOf(this.numFlexible).hashCode(),
				Boolean.valueOf(this.addWtRots).hashCode(),
				Boolean.valueOf(this.doMinimize).hashCode()
			);
		}
	}
	
	public static class ResidueFlexibility {
		
		public ArrayList<String> flexResList;
		public ArrayList<ArrayList<String>> allowedAAs;
		
		public ResidueFlexibility() {
			flexResList = new ArrayList<>();
			allowedAAs = new ArrayList<>();
		}
		
		public void addMutable(String residueNumbers, String aaNames) {
			
			// split the amino acid names
			ArrayList<String> aas = new ArrayList<>();
			for (String aaName : aaNames.split(" ")) {
				if (!aaName.isEmpty()) {
					aas.add(aaName);
				}
			}
			
			// add the residue numbers
			for (String residueNumber : residueNumbers.split(" ")) {
				if (!residueNumber.isEmpty()) {
					flexResList.add(residueNumber);
					// NOTE: for some reason, different positions can't share the same amino acid name list
					// downstream stuff just crashes for weird reasons I don't understand
					// so make sure to use a new list every time
					allowedAAs.add(new ArrayList<>(aas));
				}
			}
		}
		
		public void addFlexible(String residueNumbers) {
			addMutable(residueNumbers, "");
		}
		
		public int size() {
			return flexResList.size();
		}
	}

	public static double getRelativeError(double expected, double observed) {
		double absErr = Math.abs(expected - observed);
		if (observed == 0) {
			return absErr;
		}
		return absErr/Math.abs(observed);
	}
	
	public static Matcher<Double> isRelatively(double expected) {
		return isRelatively(expected, DefaultEpsilon);
	}
	
	public static Matcher<Double> isRelatively(final double expected, final double epsilon) {
		return new BaseMatcher<Double>() {

			@Override
			public boolean matches(Object obj) {
				double observed = ((Double)obj).doubleValue();
				if (Double.isNaN(observed)) {
					return false;
				}
				return getRelativeError(expected, observed) <= epsilon;
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
					.appendText(" that's greater than epsilon ").appendValue(epsilon);
			}
		};
	}
	
	public static double getRelativeError(BigDecimal expected, BigDecimal observed) {
		//return Math.abs(expected - observed)/Math.abs(observed);
		// NOTE: don't divide BigDecimals, since they can't represent all rationals
		BigDecimal absErr = expected.subtract(observed).abs();
		double denom = observed.abs().doubleValue();
		if (denom == 0) {
			return absErr.doubleValue();
		}
		return absErr.doubleValue()/denom;
	}
	
	public static Matcher<BigDecimal> isRelatively(BigDecimal expected) {
		return isRelatively(expected, DefaultEpsilon);
	}
	
	public static Matcher<BigDecimal> isRelatively(final BigDecimal expected, final double epsilon) {
		return new BaseMatcher<BigDecimal>() {

			@Override
			public boolean matches(Object obj) {
				BigDecimal observed = (BigDecimal)obj;
				return getRelativeError(expected, observed) <= epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				BigDecimal observed = (BigDecimal)obj;
				double relErr = getRelativeError(expected, observed);
				desc.appendText("value ").appendValue(observed.doubleValue())
					.appendText(" has relative err ").appendValue(relErr)
					.appendText(" that's greater than epsilon ").appendValue(epsilon);
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
		
		// load residue entropies
        EnvironmentVars.resTemplates.loadResEntropy("ResEntropy.dat");
	}
	
	protected static SearchProblem makeSearchProblem(EnergyMatrixConfig emConfig) {
		
		// make the search problem
		ArrayList<String> flexRes = new ArrayList<>();
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<emConfig.numFlexible; i++) {
			flexRes.add(Integer.toString(i + 1));
			allowedAAs.add(new ArrayList<String>());
		}
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;                
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", emConfig.pdbPath, 
			flexRes, allowedAAs, addWt, emConfig.doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, emConfig.addWtRots, null, 
                        false, new ArrayList<>()
		);
		
		// calculate the energy matrix, but check the cache first
		search.emat = m_energyMatrixCache.get(emConfig);
		if (search.emat == null) {
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			m_energyMatrixCache.put(emConfig, search.emat);
		}
		
		// calculate an "identity" pruning matrix (ie, no pruning)
		search.pruneMat = new PruningMatrix(search.confSpace, search.emat.getPruningInterval());
		
		return search;
	}
}
