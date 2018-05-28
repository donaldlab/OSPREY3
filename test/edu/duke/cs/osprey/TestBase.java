package edu.duke.cs.osprey;

import static org.junit.Assert.*;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.MathTools;
import org.hamcrest.BaseMatcher;
import org.hamcrest.Description;
import org.hamcrest.Matcher;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.Protractor;
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

		public void sortPositions() {
			
			class Pair {
				
				String res;
				ArrayList<String> aas;
				
				Pair(String res, ArrayList<String> aas) {
					this.res = res;
					this.aas = aas;
				}
			}
			
			List<Pair> pairs = new ArrayList<>(flexResList.size());
			for (int i=0; i<flexResList.size(); i++) {
				pairs.add(new Pair(flexResList.get(i), allowedAAs.get(i)));
			}
			Collections.sort(pairs, (Pair a, Pair b) -> {
				return Integer.parseInt(a.res) - Integer.parseInt(b.res);
			});
			
			for (int i=0; i<pairs.size(); i++) {
				Pair pair = pairs.get(i);
				flexResList.set(i, pair.res);
				allowedAAs.set(i, pair.aas);
			}	
		}
	}
	
	public static void assertConfSpacesMatch(ConfSpace confSpace, SimpleConfSpace simpleConfSpace) {
		assertThat(simpleConfSpace.positions.size(), is(confSpace.numPos));
		for (int pos=0; pos<confSpace.numPos; pos++) {
			PositionConfSpace oldpos = confSpace.posFlex.get(pos);
			Position newpos = simpleConfSpace.positions.get(pos);
			assertThat(newpos.resConfs.size(), is(oldpos.RCs.size()));
			for (int rc=0; rc<oldpos.RCs.size(); rc++) {
				RC oldrc = oldpos.RCs.get(rc);
				ResidueConf newrc = newpos.resConfs.get(rc);
				assertThat(newrc.template.name, is(oldrc.AAType));
				if (oldrc.rotNum == -1) {
					assertThat(newrc.rotamerIndex, is(nullValue()));
				} else {
					assertThat(newrc.rotamerIndex, is(oldrc.rotNum));
				}
			}
		}
	}
	
	public static double getAbsoluteError(double expected, double observed) {
		if (expected == Double.POSITIVE_INFINITY && observed == Double.POSITIVE_INFINITY) {
			return 0;
		} else if (expected == Double.NEGATIVE_INFINITY && observed == Double.NEGATIVE_INFINITY) {
			return 0;
		} else if (Double.isInfinite(expected) || Double.isInfinite(observed)) {
			return Double.POSITIVE_INFINITY;
		}
		return Math.abs(expected - observed);
	}

	public static double getRelativeError(double expected, double observed) {
		double absErr = getAbsoluteError(expected, observed);
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
	
	public static Matcher<Double> isAbsolutely(double expected) {
		return isAbsolutely(expected, DefaultEpsilon);
	}
	
	public static Matcher<Double> isAbsolutely(final double expected, final double epsilon) {
		return new BaseMatcher<Double>() {

			@Override
			public boolean matches(Object obj) {
				double observed = ((Double)obj).doubleValue();
				if (Double.isNaN(observed)) {
					return false;
				}
				return getAbsoluteError(expected, observed) <= epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				double observed = ((Double)obj).doubleValue();
				double absErr = getAbsoluteError(expected, observed);
				desc.appendText("value ").appendValue(observed)
					.appendText(" has absolute err ").appendValue(absErr)
					.appendText(" that's greater than epsilon ").appendValue(epsilon);
			}
		};
	}
	
	public static Matcher<double[]> isAbsolutely(double[] expected) {
		return isAbsolutely(expected, DefaultEpsilon);
	}
	
	public static Matcher<double[]> isAbsolutely(final double[] expected, final double epsilon) {
		return new BaseMatcher<double[]>() {
			
			int n = expected.length;

			@Override
			public boolean matches(Object obj) {
				double[] observed = (double[])obj;
				if (observed.length != n) {
					return false;
				}
				for (int i=0; i<n; i++) {
					if (getAbsoluteError(expected[i], observed[i]) > epsilon) {
						return false;
					}
				}
				return true;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				double[] observed = (double[])obj;
				
				// get the max err
				double maxAbsErr = 0;
				for (int i=0; i<n; i++) {
					maxAbsErr = Math.max(maxAbsErr, getAbsoluteError(expected[i], observed[i]));
				}
				desc.appendText("value ").appendValue(observed)
					.appendText(" has mas absolute err ").appendValue(maxAbsErr)
					.appendText(" that's greater than epsilon ").appendValue(epsilon);
			}
		};
	}
	
	public static Matcher<Double> isDegrees(double expected) {
		return isDegrees(expected, DefaultEpsilon);
	}
	
	public static Matcher<Double> isDegrees(final double expected, final double epsilon) {
		return new BaseMatcher<Double>() {
			
			private double getAbsoluteErrorDegrees(double observed) {
				// can't compare on R^1, have to compare on S^1
				return Protractor.normalizeDegrees(expected - observed);
			}

			@Override
			public boolean matches(Object obj) {
				double observed = ((Double)obj).doubleValue();
				if (Double.isNaN(observed)) {
					return false;
				}
				return getAbsoluteErrorDegrees(observed) <= epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				double observed = ((Double)obj).doubleValue();
				double absErr = getAbsoluteErrorDegrees(observed);
				desc.appendText("value ").appendValue(observed)
					.appendText(" has absolute err ").appendValue(absErr)
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
	
	
	@Deprecated
	protected static ForcefieldParams makeDefaultFFParams() {
		return new ForcefieldParams();
	}
	
	@Deprecated
	protected static void initDefaultEnvironment() {
		
		// make energy function
		EnvironmentVars.curEFcnGenerator = new EnergyFunctionGenerator(new ForcefieldParams());
		
		// make residue templates
		EnvironmentVars.resTemplates = new ResidueTemplateLibrary.Builder().build();
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

	public static void fileForWriting(String filename, Consumer<File> block) {
		File file = new File(filename);
		file.delete();
		try {
			block.accept(file);
		} finally {
			file.delete();
		}
	}
}
