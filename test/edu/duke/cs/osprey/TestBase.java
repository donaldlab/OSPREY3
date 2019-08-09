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
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.BigExp;
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
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
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

				if (observed.length != n) {

					desc.appendText("value ").appendValue(observed)
						.appendText(" has length ").appendValue(observed.length);

				} else {
				
					// get the max err
					double maxAbsErr = 0;
					for (int i=0; i<n; i++) {
						maxAbsErr = Math.max(maxAbsErr, getAbsoluteError(expected[i], observed[i]));
					}
					desc.appendText("value ").appendValue(observed)
						.appendText(" has max absolute err ").appendValue(maxAbsErr)
						.appendText(" that's greater than epsilon ").appendValue(epsilon);
				}
			}
		};
	}


	private abstract static class EpsilonApplier<T> {

		public final double epsilon;

		protected EpsilonApplier(double epsilon) {
			this.epsilon = epsilon;
		}

		abstract T apply(T thing);
		abstract String term();

		static EpsilonApplier<DoubleBounds> doubleBoundsAbsolute(double epsilon) {
			return new EpsilonApplier<DoubleBounds>(epsilon) {
				@Override
				public String term() {
					return "absolutely";
				}
				@Override
				public DoubleBounds apply(DoubleBounds bounds) {
					return new DoubleBounds(
						bounds.lower - epsilon,
						bounds.upper + epsilon
					);
				}
			};
		}

		static EpsilonApplier<DoubleBounds> doubleBoundsRelative(double epsilon) {
			return new EpsilonApplier<DoubleBounds>(epsilon) {
				@Override
				public String term() {
					return "relatively";
				}
				@Override
				public DoubleBounds apply(DoubleBounds bounds) {
					double loFactor = 1.0 + (bounds.lower > 0 ? -epsilon : +epsilon);
					double hiFactor = 1.0 + (bounds.upper > 0 ? +epsilon : -epsilon);
					return new DoubleBounds(
						bounds.lower*loFactor,
						bounds.upper*hiFactor
					);
				}
			};
		}

		static EpsilonApplier<BigDecimalBounds> bigDecimalBoundsAbsolute(double epsilon) {
			return new EpsilonApplier<BigDecimalBounds>(epsilon) {

				// use a math context with fixed precision to keep additions from being really slow!
				MathContext mathContext = new MathContext(32, RoundingMode.HALF_UP);
				BigDecimal bigEpsilon = MathTools.biggen(epsilon);

				@Override
				public String term() {
					return "absolutely";
				}
				@Override
				public BigDecimalBounds apply(BigDecimalBounds bounds) {
					return new BigDecimalBounds(
						bounds.lower.subtract(bigEpsilon, mathContext),
						bounds.upper.add(bigEpsilon, mathContext)
					);
				}
			};
		}

		static EpsilonApplier<BigDecimalBounds> bigDecimalBoundsRelative(double epsilon) {
			return new EpsilonApplier<BigDecimalBounds>(epsilon) {

				@Override
				public String term() {
					return "relatively";
				}
				@Override
				public BigDecimalBounds apply(BigDecimalBounds bounds) {
					// multiplications are always fast enough, even without a math context
					BigDecimal loFactor = MathTools.biggen(1.0 + (MathTools.isPositive(bounds.lower) ? -epsilon : +epsilon));
					BigDecimal hiFactor = MathTools.biggen(1.0 + (MathTools.isPositive(bounds.upper) ? +epsilon : -epsilon));
					return new BigDecimalBounds(
						bounds.lower.multiply(loFactor),
						bounds.upper.multiply(hiFactor)
					);
				}
			};
		}

		static EpsilonApplier<BigExp.Bounds> bigExpBoundsAbsolute(double epsilon) {
			return new EpsilonApplier<BigExp.Bounds>(epsilon) {

				BigExp bigEpsilon = new BigExp(epsilon);

				@Override
				public String term() {
					return "absolutely";
				}
				@Override
				public BigExp.Bounds apply(BigExp.Bounds bounds) {
					bounds = new BigExp.Bounds(new BigExp(bounds.lower), new BigExp(bounds.upper));
					bounds.lower.sub(bigEpsilon);
					bounds.upper.add(bigEpsilon);
					return bounds;
				}
			};
		}

		static EpsilonApplier<BigExp.Bounds> bigExpBoundsRelative(double epsilon) {
			return new EpsilonApplier<BigExp.Bounds>(epsilon) {

				@Override
				public String term() {
					return "relatively";
				}
				@Override
				public BigExp.Bounds apply(BigExp.Bounds bounds) {
					// multiplications are always fast enough, even without a math context
					double loFactor = 1.0 + (bounds.lower.isPositive() ? -epsilon : +epsilon);
					double hiFactor = 1.0 + (bounds.upper.isPositive() ? +epsilon : -epsilon);
					bounds = new BigExp.Bounds(new BigExp(bounds.lower), new BigExp(bounds.upper));
					bounds.lower.mult(loFactor);
					bounds.upper.mult(hiFactor);
					return bounds;
				}
			};
		}
	}

	private static class BigDecimalBoundsMatchingBigDecimal extends BaseMatcher<BigDecimalBounds> {

		private EpsilonApplier<BigDecimalBounds> epsilonApplier;
		private BigDecimal expected;

		public BigDecimalBoundsMatchingBigDecimal(EpsilonApplier<BigDecimalBounds> epsilonApplier, BigDecimal expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			BigDecimalBounds observed = (BigDecimalBounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			BigDecimalBounds observed = (BigDecimalBounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<BigDecimalBounds> isAbsoluteBound(BigDecimal expected, double epsilon) {
		return new BigDecimalBoundsMatchingBigDecimal(EpsilonApplier.bigDecimalBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<BigDecimalBounds> isRelativeBound(BigDecimal expected, double epsilon) {
		return new BigDecimalBoundsMatchingBigDecimal(EpsilonApplier.bigDecimalBoundsRelative(epsilon), expected);
	}

	private static class BigDecimalBoundsMatchingBigDecimalBounds extends BaseMatcher<BigDecimalBounds> {

		EpsilonApplier<BigDecimalBounds> epsilonApplier;
		BigDecimalBounds expected;

		public BigDecimalBoundsMatchingBigDecimalBounds(EpsilonApplier<BigDecimalBounds> epsilonApplier, BigDecimalBounds expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			BigDecimalBounds observed = (BigDecimalBounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			BigDecimalBounds observed = (BigDecimalBounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<BigDecimalBounds> isAbsoluteBound(BigDecimalBounds expected, double epsilon) {
		return new BigDecimalBoundsMatchingBigDecimalBounds(EpsilonApplier.bigDecimalBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<BigDecimalBounds> isRelativeBound(BigDecimalBounds expected, double epsilon) {
		return new BigDecimalBoundsMatchingBigDecimalBounds(EpsilonApplier.bigDecimalBoundsRelative(epsilon), expected);
	}

	private static class BigExpBoundsMatchingBigExp extends BaseMatcher<BigExp.Bounds> {

		private EpsilonApplier<BigExp.Bounds> epsilonApplier;
		private BigExp expected;

		public BigExpBoundsMatchingBigExp(EpsilonApplier<BigExp.Bounds> epsilonApplier, BigExp expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			BigExp.Bounds observed = (BigExp.Bounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			BigExp.Bounds observed = (BigExp.Bounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<BigExp.Bounds> isAbsoluteBound(BigExp expected, double epsilon) {
		return new BigExpBoundsMatchingBigExp(EpsilonApplier.bigExpBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<BigExp.Bounds> isRelativeBound(BigExp expected, double epsilon) {
		return new BigExpBoundsMatchingBigExp(EpsilonApplier.bigExpBoundsRelative(epsilon), expected);
	}

	private static class BigExpBoundsMatchingBigExpBounds extends BaseMatcher<BigExp.Bounds> {

		EpsilonApplier<BigExp.Bounds> epsilonApplier;
		BigExp.Bounds expected;

		public BigExpBoundsMatchingBigExpBounds(EpsilonApplier<BigExp.Bounds> epsilonApplier, BigExp.Bounds expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			BigExp.Bounds observed = (BigExp.Bounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			BigExp.Bounds observed = (BigExp.Bounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<BigExp.Bounds> isAbsoluteBound(BigExp.Bounds expected, double epsilon) {
		return new BigExpBoundsMatchingBigExpBounds(EpsilonApplier.bigExpBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<BigExp.Bounds> isRelativeBound(BigExp.Bounds expected, double epsilon) {
		return new BigExpBoundsMatchingBigExpBounds(EpsilonApplier.bigExpBoundsRelative(epsilon), expected);
	}

	private static class DoubleBoundsMatchingDouble extends BaseMatcher<DoubleBounds> {

		private EpsilonApplier<DoubleBounds> epsilonApplier;
		private double expected;

		public DoubleBoundsMatchingDouble(EpsilonApplier<DoubleBounds> epsilonApplier, double expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			DoubleBounds observed = (DoubleBounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			DoubleBounds observed = (DoubleBounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<DoubleBounds> isAbsoluteBound(double expected, double epsilon) {
		return new DoubleBoundsMatchingDouble(EpsilonApplier.doubleBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<DoubleBounds> isRelativeBound(double expected, double epsilon) {
		return new DoubleBoundsMatchingDouble(EpsilonApplier.doubleBoundsRelative(epsilon), expected);
	}

	private static class DoubleBoundsMatchingDoubleBounds extends BaseMatcher<DoubleBounds> {

		private EpsilonApplier<DoubleBounds> epsilonApplier;
		private DoubleBounds expected;

		public DoubleBoundsMatchingDoubleBounds(EpsilonApplier<DoubleBounds> epsilonApplier, DoubleBounds expected) {
			this.epsilonApplier = epsilonApplier;
			this.expected = expected;
		}

		@Override
		public boolean matches(Object obj) {
			DoubleBounds observed = (DoubleBounds)obj;
			return observed.isValid() && epsilonApplier.apply(observed).contains(expected);
		}

		@Override
		public void describeTo(Description desc) {
			desc.appendText("bounds ").appendValue(expected)
				.appendText(" " + epsilonApplier.term() + " within epsilon ").appendValue(epsilonApplier.epsilon);
		}

		@Override
		public void describeMismatch(Object obj, Description desc) {
			DoubleBounds observed = (DoubleBounds)obj;
			if (!observed.isValid()) {
				desc.appendValue(observed).appendText(" is not a valid bound");
			} else {
				if (!epsilonApplier.apply(observed).contains(expected)) {
					desc.appendValue(observed)
						.appendText(" (with epsilon: ")
						.appendValue(epsilonApplier.apply(observed))
						.appendText(") does not bound " + epsilonApplier.term() + " within epsilon");
				}
			}
		}
	}
	public static Matcher<DoubleBounds> isAbsoluteBound(DoubleBounds expected, double epsilon) {
		return new DoubleBoundsMatchingDoubleBounds(EpsilonApplier.doubleBoundsAbsolute(epsilon), expected);
	}
	public static Matcher<DoubleBounds> isRelativeBound(DoubleBounds expected, double epsilon) {
		return new DoubleBoundsMatchingDoubleBounds(EpsilonApplier.doubleBoundsRelative(epsilon), expected);
	}


	public static Matcher<Double> isDegrees(double expected) {
		return isDegrees(expected, DefaultEpsilon);
	}
	
	public static Matcher<Double> isDegrees(final double expected, final double epsilon) {
		return new BaseMatcher<Double>() {
			
			private double getAbsoluteErrorDegrees(double observed) {
				// can't compare on R^1, have to compare on S^1
				return Protractor.getDistDegrees(expected, observed);
			}

			@Override
			public boolean matches(Object obj) {
				double observed = (Double)obj;
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
				double observed = (Double)obj;
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

	public static double getRelativeError(BigExp expected, BigExp observed) {
		BigExp err = new BigExp(expected);
		err.sub(observed);
		err.abs();
		err.div(observed);
		return err.toDouble();
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
				desc.appendText("close to ").appendValue(String.format("%e", expected));
			}
			
			@Override
			public void describeMismatch(Object obj, Description desc) {
				BigDecimal observed = (BigDecimal)obj;
				double relErr = getRelativeError(expected, observed);
				desc.appendText("value ").appendValue(String.format("%e", observed))
					.appendText(" has relative err ").appendValue(relErr)
					.appendText(" that's greater than epsilon ").appendValue(epsilon);
			}
		};
	}

	public static Matcher<BigExp> isRelatively(BigExp expected) {
		return isRelatively(expected, DefaultEpsilon);
	}

	public static Matcher<BigExp> isRelatively(final BigExp expected, final double epsilon) {
		return new BaseMatcher<BigExp>() {

			@Override
			public boolean matches(Object obj) {
				BigExp observed = (BigExp)obj;
				return getRelativeError(expected, observed) <= epsilon;
			}

			@Override
			public void describeTo(Description desc) {
				desc.appendText("close to ").appendValue(expected);
			}

			@Override
			public void describeMismatch(Object obj, Description desc) {
				BigExp observed = (BigExp)obj;
				double relErr = getRelativeError(expected, observed);
				desc.appendText("value ").appendValue(observed)
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

	/**
	 * a File subclass that automatically deletes itself when finished
	 */
	public static class TempFile extends File implements AutoCloseable {

		public TempFile(String filename) {
			super(filename);
		}

		public TempFile(File dir, String filename) {
			super(dir, filename);
		}

		@Override
		public void close() {
			delete();
		}
	}
}
