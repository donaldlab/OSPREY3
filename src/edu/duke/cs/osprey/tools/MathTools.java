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

package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;


public class MathTools {
	
	public static int divUp(int num, int denom) {
		return (num + denom - 1)/denom;
	}
	
	public static int roundUpToMultiple(int val, int base) {
		int mod = val % base;
		if (mod == 0) {
			return val;
		}
		return val + base - mod;
	}

	public static <T> List<List<T>> powerset(List<T> list) {
		return powersetUpTo(list, list.size());
	}

	public static <T> List<List<T>> powersetUpTo(List<T> list, int size) {

		// adapted from powerset algorithm:
		// http://rosettacode.org/wiki/Power_set#Iterative

		List<List<T>> powerset = new ArrayList<>();

		// start with the empty set
		powerset.add(new ArrayList<>());

		for (T item : list) {

			for (List<T> subset : new ArrayList<>(powerset)) {

				// expand subset if within size
				if (subset.size() < size) {
					List<T> newSubset = new ArrayList<>(subset);
					newSubset.add(item);
					powerset.add(newSubset);
				}
			}
		}
		return powerset;
	}

	public static <T> Iterable<List<T>> cartesianProduct(List<List<T>> lists) {

		// stack overflow is awesome for the lazy programmer:
		// https://stackoverflow.com/questions/714108/cartesian-product-of-arbitrary-sets-in-java
		// I mean, sure I could figure out how to do this stuff from first principles, but why would I want to?

		return () -> {

			return new Iterator<List<T>>() {

				private final int n;
				private boolean hasNext;
				private final int[] lengths;
				private final int[] indices;

				{
					n = lists.size();
					hasNext = n > 0;

					lengths = new int[n];
					for (int i=0; i<n; i++) {
						lengths[i] = lists.get(i).size();
						if (lengths[i] <= 0) {
							hasNext = false;
						}
					}

					indices = new int[n];
					Arrays.fill(indices, 0);
				}

				public boolean hasNext() {
					return hasNext;
				}

				public List<T> next() {

					// build the output
					List<T> result = new ArrayList<>();
					for (int i=0; i<n; i++) {
						result.add(lists.get(i).get(indices[i]));
					}

					// advance the indices
					for (int i = n - 1; i >= 0; i--) {
						if (indices[i] == lengths[i] - 1) {
							indices[i] = 0;
							if (i == 0) {
								hasNext = false;
							}
						} else {
							indices[i]++;
							break;
						}
					}

					return result;
				}
			};
		};
	}

	public static boolean isZero(BigInteger i) {
		return i.compareTo(BigInteger.ZERO) == 0;
	}

	// HACKHACK: need a way to represent infinity for big decimals
	// our "infinity" value still needs an actual value though,
	// so use the biggest number we can and hope it never happens in a real design

	// TODO: we should really get rid of these hacks
	// and use a better arbitrary-precision float lib instead of BigDecimal
	// maybe something like JScience would be better:
	// http://jscience.org/api/org/jscience/mathematics/number/package-summary.html

	public static class MagicBigDecimal extends BigDecimal {

		public final double doubleValue;

		public MagicBigDecimal(double doubleValue) {

			// pick an arbitrary value in the range of BigDecimal values
			// that we can reserve for our magic numbers
			super(new BigInteger("1"), Integer.MAX_VALUE);

			this.doubleValue = doubleValue;
		}

		@Override
		public double doubleValue() {
			return doubleValue;
		}

		@Override
		public String toString() {
			return Double.toString(doubleValue);
		}

		// try to enforce the prohibition against doing math with our magic values
		// these checks won't cover all our bases, but it's a start
		@Override public BigDecimal add(BigDecimal other) { throw new DontDoRawMathWithMagicException(); }
		@Override public BigDecimal add(BigDecimal other, MathContext context) { throw new DontDoRawMathWithMagicException(); }

		@Override public BigDecimal subtract(BigDecimal other) { throw new DontDoRawMathWithMagicException(); }
		@Override public BigDecimal subtract(BigDecimal other, MathContext context) { throw new DontDoRawMathWithMagicException(); }

		@Override public BigDecimal multiply(BigDecimal other) { throw new DontDoRawMathWithMagicException(); }
		@Override public BigDecimal multiply(BigDecimal other, MathContext context) { throw new DontDoRawMathWithMagicException(); }

		@Override public BigDecimal divide(BigDecimal other) { throw new DontDoRawMathWithMagicException(); }
		@Override public BigDecimal divide(BigDecimal other, MathContext context) { throw new DontDoRawMathWithMagicException(); }
		@Override public BigDecimal divide(BigDecimal divisor, int scale, int roundingMode) { throw new DontDoRawMathWithMagicException(); }

		@Override public int precision() { throw new DontDoRawMathWithMagicException(); }
		@Override public int scale() { throw new DontDoRawMathWithMagicException(); }
		@Override public int signum() { throw new DontDoRawMathWithMagicException(); }

		public static class DontDoRawMathWithMagicException extends UnsupportedOperationException {}

		@Override
		public boolean equals(Object other) {
			// these are only singleton instances
			return this == other;
		}
	}

	/**
	 * Sadly, Java's BigDecimal can't encode values of +Infinity or -Infinity =(
	 * These magic constants are a complete hack to try to work around that.
	 * Don't try to use this value in any arithmetic. That won't work.
	 * Just compare BigDecimal references to check for infinity, e.g. if (myval == BigPositiveInfinity) { ... }
	 */
	public static final BigDecimal BigPositiveInfinity = new MagicBigDecimal(Double.POSITIVE_INFINITY);

	/** See BigPositiveInfinity for usage instructions */
	public static final BigDecimal BigNegativeInfinity = new MagicBigDecimal(Double.NEGATIVE_INFINITY);

	public static final BigDecimal BigNaN = new MagicBigDecimal(Double.NaN);

	public static BigDecimal biggen(double val) {
		if (val == Double.POSITIVE_INFINITY) {
			return BigPositiveInfinity;
		} else if (val == Double.NEGATIVE_INFINITY) {
			return BigNegativeInfinity;
		} else if (Double.isNaN(val)) {
			return BigNaN;
		} else if (val == 0.0) {
			return BigDecimal.ZERO;
		} else if (val == 1.0) {
			return BigDecimal.ONE;
		} else {
			return BigDecimal.valueOf(val);
		}
	}

	public static BigDecimal biggen(long val) {
		return new BigDecimal(val);
	}

	public static BigDecimalBounds biggen(double lower, double upper) {
		return new BigDecimalBounds(biggen(lower), biggen(upper));
	}

	public static int compare(BigDecimal a, BigDecimal b) {
		// a < b => -1
		// a == b => 0
		// a > b => 1
		if (a == BigNaN || b == BigNaN) {
			throw new IllegalArgumentException("can't compare NaN");
		}
		if (a == BigPositiveInfinity) {
			if (b == BigPositiveInfinity) {
				return 0;
			} else {
				return 1;
			}
		} else if (a == BigNegativeInfinity) {
			if (b == BigNegativeInfinity) {
				return 0;
			} else {
				return -1;
			}
		} else {
			if (b == BigPositiveInfinity) {
				return -1;
			} else if (b == BigNegativeInfinity) {
				return 1;
			} else {
				return a.compareTo(b);
			}
		}
	}

	/**
	 * Tests for sameness of values,
	 * rather than BigDecimal.equals(), which tests sameness of representation
	 * (i.e., the same value can have multiple representations in BigDecimal)
	 */
	public static boolean isSameValue(BigDecimal a, BigDecimal b) {
		return a == b || compare(a, b) == 0;
	}

	public static boolean isAbsolutelySame(BigDecimal a, BigDecimal b, double epsilon) {
		return a == b || a.subtract(b).abs().doubleValue() <= epsilon;
	}

	public static boolean isRelativelySame(BigDecimal a, BigDecimal b, MathContext context, double epsilon) {
		if (a == b) {
			return true;
		} else if (a instanceof MagicBigDecimal || b instanceof MagicBigDecimal) {
			return false;
		}
		BigDecimal numerator = a.subtract(b).abs();
		BigDecimal denominator = a.abs();
		if (isZero(numerator) && isZero(denominator)) {
			return true;
		} else if (isZero(numerator)) {
			return true;
		} else if (isZero(denominator)) {
			return false;
		} else {
			return numerator.divide(denominator, context).doubleValue() <= epsilon;
		}
	}

	public static boolean isZero(BigDecimal d) {
		return d == BigDecimal.ZERO || (isFinite(d) && d.compareTo(BigDecimal.ZERO) == 0);
	}

	public static boolean isPositive(BigDecimal d) {
		if (d == BigPositiveInfinity) {
			return true;
		} else if (d == BigNegativeInfinity) {
			return false;
		} else if (d == BigNaN) {
			return false;
		} else {
			return isGreaterThan(d, BigDecimal.ZERO);
		}
	}

	public static boolean isNegative(BigDecimal d) {
		if (d == BigPositiveInfinity) {
			return false;
		} else if (d == BigNegativeInfinity) {
			return true;
		} else if (d == BigNaN) {
			return false;
		} else {
			return isLessThan(d, BigDecimal.ZERO);
		}
	}

	public static boolean isInf(BigDecimal d) {
		return d == BigPositiveInfinity || d == BigNegativeInfinity;
	}

	public static boolean isFinite(BigDecimal d) {
		return !isInf(d) && !isNaN(d);
	}

	public static boolean isNaN(BigDecimal d) {
		return d == BigNaN;
	}

	/** return a < b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isLessThan(BigDecimal a, BigDecimal b) {
		return compare(a, b) == -1;
	}

	/** return a <= b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isLessThanOrEqual(BigDecimal a, BigDecimal b) {
		return isLessThan(a, b) || isSameValue(a, b);
	}

	/** return a > b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isGreaterThan(BigDecimal a, BigDecimal b) {
		return compare(a, b) == 1;
	}

	/** return a >= b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isGreaterThanOrEqual(BigDecimal a, BigDecimal b) {
		return isGreaterThan(a, b) || isSameValue(a, b);
	}

	/** return -a, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigNegate(BigDecimal a) {
		if (a == BigNaN) {
			return BigNaN;
		} else if (a == BigNegativeInfinity) {
			return BigPositiveInfinity;
		} else if (a == BigPositiveInfinity) {
			return BigNegativeInfinity;
		} else {
			return a.negate();
		}
	}

	/** return a + b, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigAdd(BigDecimal a, BigDecimal b, MathContext context) {
		if (a == BigNaN || b == BigNaN) {
			return BigNaN;
		}
		if (a == BigPositiveInfinity) {
			if (b == BigNegativeInfinity) {
				return BigNaN;
			} else {
				return BigPositiveInfinity;
			}
		} else if (a == BigNegativeInfinity) {
			if (b == BigPositiveInfinity) {
				return BigNaN;
			} else {
				return BigNegativeInfinity;
			}
		} else {
			if (b == BigPositiveInfinity) {
				return BigPositiveInfinity;
			} else if (b == BigNegativeInfinity) {
				return BigNegativeInfinity;
			} else {
				return a.add(b, context);
			}
		}
	}

	/** return a - b, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigSubtract(BigDecimal a, BigDecimal b, MathContext context) {
		if (a == BigNaN || b == BigNaN) {
			return BigNaN;
		}
		if (a == BigPositiveInfinity) {
			if (b == BigPositiveInfinity) {
				return BigNaN;
			} else {
				return BigPositiveInfinity;
			}
		} else if (a == BigNegativeInfinity) {
			if (b == BigNegativeInfinity) {
				return BigNaN;
			} else {
				return BigNegativeInfinity;
			}
		} else {
			if (b == BigPositiveInfinity) {
				return BigNegativeInfinity;
			} else if (b == BigNegativeInfinity) {
				return BigPositiveInfinity;
			} else {
				return a.subtract(b, context);
			}
		}
	}

	/** return a*b, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigMultiply(BigDecimal a, BigDecimal b, MathContext context) {
		if (a == BigNaN) {
			return BigNaN;
		} else if (a == BigPositiveInfinity) {
			if (b == BigNaN || b == BigNegativeInfinity) {
				return BigNaN;
			} else if (b == BigPositiveInfinity) {
				return BigPositiveInfinity;
			} else {
				int cmp = b.compareTo(BigDecimal.ZERO);
				if (cmp == 0) {
					return BigDecimal.ZERO;
				} else if (cmp > 0) {
					return BigPositiveInfinity;
				} else {
					return BigNegativeInfinity;
				}
			}
		} else if (a == BigNegativeInfinity) {
			if (b == BigNaN || b == BigPositiveInfinity) {
				return BigNaN;
			} else if (b == BigNegativeInfinity) {
				return BigPositiveInfinity;
			} else {
				int cmp = b.compareTo(BigDecimal.ZERO);
				if (cmp == 0) {
					return BigDecimal.ZERO;
				} else if (cmp > 0) {
					return BigNegativeInfinity;
				} else {
					return BigPositiveInfinity;
				}
			}
		} else {
			if (b == BigNaN) {
				return BigNaN;
			} else if (b == BigPositiveInfinity) {
				int cmp = a.compareTo(BigDecimal.ZERO);
				if (cmp == 0) {
					return BigDecimal.ZERO;
				} else if (cmp > 0) {
					return BigPositiveInfinity;
				} else {
					return BigNegativeInfinity;
				}
			} else if (b == BigNegativeInfinity) {
				int cmp = a.compareTo(BigDecimal.ZERO);
				if (cmp == 0) {
					return BigDecimal.ZERO;
				} else if (cmp > 0) {
					return BigNegativeInfinity;
				} else {
					return BigPositiveInfinity;
				}
			} else {
				try {
					return a.multiply(b, context);
				} catch (ArithmeticException ex) {
					if (ex.getMessage().equals("Underflow")) {
						//log("WARN: multiplication underflows. Assuming %.4e * %.4e is just 0", a, b);
						return BigDecimal.ZERO;
					} else {
						throw ex;
					}
				}
			}
		}
	}

	/** return a/b, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigDivide(BigDecimal a, BigDecimal b, MathContext context) {
		if (a == BigNaN || b == BigNaN) {
			return BigNaN;
		}
		if (a == BigPositiveInfinity) {
			if (isInf(b)) {
				return BigNaN;
			} else if (b.signum() >= 0) {
				return BigPositiveInfinity;
			} else {
				return BigNegativeInfinity;
			}
		} else if (a == BigNegativeInfinity) {
			if (isInf(b)) {
				return BigNaN;
			} else if (b.signum() >= 0) {
				return BigNegativeInfinity;
			} else {
				return BigPositiveInfinity;
			}
		} else {
			if (isInf(b)) {
				return BigDecimal.ZERO;
			} else if (isZero(a) && isZero(b)) {
				return BigNaN;
			} else if (isZero(b)) {
				return BigPositiveInfinity;
			} else {
				return a.divide(b, context);
			}
		}
	}

	/** return a/b/c, correctly handling -Inf, +Inf, and NaN */
	public static BigDecimal bigDivideDivide(BigDecimal a, BigDecimal b, BigDecimal c, MathContext context) {
		return bigDivide(bigDivide(a, b, context), c, context);
	}

	/**
	 * returns log10(x+1)
	 *
	 * chosen because it maps [0,inf] to [0,inf]
	 * whereas just log10 maps [0,inf] to [-inf,inf]
	 **/
	public static double log10p1(BigDecimal x) {
		if (x == BigPositiveInfinity) {
			return Double.POSITIVE_INFINITY;
		} else if (x == BigNegativeInfinity || x == BigNaN) {
			return Double.NaN;
		} else {
			return Math.log10(x.add(BigDecimal.ONE).doubleValue());
		}
	}

	public static double log10p1(double x) {
		return Math.log10(x + 1);
	}

	public static String formatBytes(long bytes) {
		if (bytes < 1024) {
			return String.format("%d B", bytes);
		}
		double kibibytes = (double)bytes/1024;
		if (kibibytes < 128) {
			return String.format("%.1f KiB", kibibytes);
		} else if (kibibytes < 1024) {
			return String.format("%.0f KiB", kibibytes);
		}
		double mebibytes = kibibytes/1024;
		if (mebibytes < 128) {
			return String.format("%.1f MiB", mebibytes);
		} else if (mebibytes < 1024) {
			return String.format("%.0f MiB", mebibytes);
		}
		double gibibytes = mebibytes/1024;
		if (gibibytes < 128) {
			return String.format("%.1f GiB", gibibytes);
		} else if (gibibytes < 1024) {
			return String.format("%.0f GiB", gibibytes);
		}
		double tebibytes = gibibytes/1024;
		if (tebibytes < 128) {
			return String.format("%.1f TiB", tebibytes);
		} else if (tebibytes < 1024) {
			return String.format("%.0f TiB", tebibytes);
		}
		double pebibytes = tebibytes/1024;
		if (pebibytes < 128) {
			return String.format("%.1f PiB", pebibytes);
		} else {
			return String.format("%.0f PiB", pebibytes);
		}
	}

	public static enum Optimizer {

		Minimize {

			@Override
			public float initFloat() {
				return Float.POSITIVE_INFINITY;
			}

			@Override
			public double initDouble() {
				return Double.POSITIVE_INFINITY;
			}

			@Override
			public int initInt() {
				return Integer.MAX_VALUE;
			}

			@Override
			public long initLong() {
				return Long.MAX_VALUE;
			}

			@Override
			public BigDecimal initBigDecimal() {
				return MathTools.BigPositiveInfinity;
			}

			@Override
			public boolean isBetter(float newval, float oldval) {
				return newval < oldval;
			}

			@Override
			public boolean isBetter(double newval, double oldval) {
				return newval < oldval;
			}

			@Override
			public boolean isBetter(int newval, int oldval) {
				return newval < oldval;
			}

			@Override
			public boolean isBetter(long newval, long oldval) {
				return newval < oldval;
			}

			@Override
			public boolean isBetter(BigDecimal newval, BigDecimal oldval) {
				return MathTools.isLessThan(newval, oldval);
			}

			@Override
			public Optimizer reverse() {
				return Maximize;
			}
		},

		Maximize {

			@Override
			public float initFloat() {
				return Float.NEGATIVE_INFINITY;
			}

			@Override
			public double initDouble() {
				return Double.NEGATIVE_INFINITY;
			}

			@Override
			public int initInt() {
				return Integer.MIN_VALUE;
			}

			@Override
			public long initLong() {
				return Long.MIN_VALUE;
			}

			@Override
			public BigDecimal initBigDecimal() {
				return MathTools.BigNegativeInfinity;
			}

			@Override
			public boolean isBetter(float newval, float oldval) {
				return newval > oldval;
			}

			@Override
			public boolean isBetter(double newval, double oldval) {
				return newval > oldval;
			}

			@Override
			public boolean isBetter(int newval, int oldval) {
				return newval > oldval;
			}

			@Override
			public boolean isBetter(long newval, long oldval) {
				return newval > oldval;
			}

			@Override
			public boolean isBetter(BigDecimal newval, BigDecimal oldval) {
				return MathTools.isGreaterThan(newval, oldval);
			}

			@Override
			public Optimizer reverse() {
				return Minimize;
			}
		};

		public abstract float initFloat();
		public abstract double initDouble();
		public abstract int initInt();
		public abstract long initLong();
		public abstract BigDecimal initBigDecimal();

		public float opt(float a, float b) {
			return isBetter(a, b) ? a : b;
		}

		public double opt(double a, double b) {
			return isBetter(a, b) ? a : b;
		}

		public int opt(int a, int b) {
			return isBetter(a, b) ? a : b;
		}

		public long opt(long a, long b) {
			return isBetter(a, b) ? a : b;
		}

		public BigDecimal opt(BigDecimal a, BigDecimal b) {
			return isBetter(a, b) ? a : b;
		}

		public abstract boolean isBetter(float newval, float oldval);
		public abstract boolean isBetter(double newval, double oldval);
		public abstract boolean isBetter(int newval, int oldval);
		public abstract boolean isBetter(long newval, long oldval);
		public abstract boolean isBetter(BigDecimal newval, BigDecimal oldval);

		public abstract Optimizer reverse();
	}


	public static class DoubleBounds {

		public double lower;
		public double upper;

		public DoubleBounds() {
			this(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
		}

		public DoubleBounds(double lower, double upper) {
			this.lower = lower;
			this.upper = upper;
		}

		public void expand(double p) {
			lower = Math.min(lower, p);
			upper = Math.max(upper, p);
		}

		public double size() {
			return upper - lower;
		}

		public boolean isValid() {
			return lower <= upper;
		}

		public boolean contains(double val) {
			return val >= lower && val <= upper;
		}

		public boolean contains(DoubleBounds val) {
			return val.lower >= this.lower && val.upper <= this.upper;
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				Double.hashCode(lower),
				Double.hashCode(upper)
			);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof DoubleBounds && equals((DoubleBounds)other);
		}

		public boolean equals(DoubleBounds other) {
			return this.lower == other.lower
				&& this.upper == other.upper;
		}

		@Override
		public String toString() {
			return toString(null, null);
		}

		public String toString(int precision) {
			return toString(precision, null);
		}

		public String toString(Integer precision, Integer width) {
			String spec = "%" + (width != null ? width : "") + (precision != null ? "." + precision : "") + "f";
			return String.format("[" + spec + "," + spec + "]", lower, upper);
		}
	}

	public static class BigDecimalBounds {

		public BigDecimal lower;
		public BigDecimal upper;

		public BigDecimalBounds() {
			this(BigNegativeInfinity, BigPositiveInfinity);
		}

		public BigDecimalBounds(BigDecimalBounds other) {
			this(other.lower, other.upper);
		}

		public BigDecimalBounds(BigDecimal lower, BigDecimal upper) {
			this.lower = lower;
			this.upper = upper;
		}

		public BigDecimalBounds(BigDecimal val) {
			this(val, val);
		}

		public BigDecimalBounds(double lower, double upper) {
			this(biggen(lower), biggen(upper));
		}

		public BigDecimal size(MathContext mathContext) {
			return new BigMath(mathContext)
				.set(upper)
				.sub(lower)
				.get();
		}

		public boolean isValid() {
			return MathTools.isGreaterThanOrEqual(upper, lower);
		}

		public boolean contains(BigDecimal d) {
			return MathTools.isGreaterThanOrEqual(d, lower)
				&& MathTools.isLessThanOrEqual(d, upper);
		}

		public boolean contains(BigDecimalBounds d) {
			return MathTools.isGreaterThanOrEqual(d.lower, lower)
				&& MathTools.isLessThanOrEqual(d.upper, upper);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof BigDecimalBounds && equals((BigDecimalBounds)other);
		}

		public boolean equals(BigDecimalBounds other) {
			return MathTools.isSameValue(this.lower, other.lower)
				&& MathTools.isSameValue(this.upper, other.upper);
		}

		@Override
		public int hashCode() {
			return HashCalculator.combineHashes(
				lower.hashCode(),
				upper.hashCode()
			);
		}

		@Override
		public String toString() {
			return String.format("[%e,%e]", lower, upper);
		}
	}

	public static class BigIntegerBounds {

		public BigInteger lower;
		public BigInteger upper;

		public BigIntegerBounds(BigInteger lower, BigInteger upper) {
			this.lower = lower;
			this.upper = upper;
		}

		@Override
		public String toString() {
			return String.format("[%d,%d]", lower, upper);
		}
	}

	public static class GridIterable implements Iterable<int[]> {

		public final int[] dimensions;

		public GridIterable(int[] dimensions) {

			for (int d : dimensions) {
				if (d <= 0) {
					throw new IllegalArgumentException("invalid dimensions: " + Arrays.toString(dimensions));
				}
			}

			this.dimensions = dimensions;
		}

		@Override
		public Iterator<int[]> iterator() {
			return new Iterator<int[]>() {

				int[] indices = new int[dimensions.length];
				boolean hasNext = true;

				{
					// start indices at one pos before all zeros,
					// so the first call to next() moves to all zeros
					Arrays.fill(indices, 0);
					indices[0] = -1;
				}

				@Override
				public boolean hasNext() {
					return hasNext;
				}

				@Override
				public int[] next() {

					if (!hasNext) {
						throw new NoSuchElementException();
					}

					// advance to the next index
					indices[0]++;
					for (int d=0; d<dimensions.length; d++) {
						if (indices[d] >= dimensions[d]) {
							if (d + 1 == dimensions.length) {
								throw new UnpossibleError();
							}
							indices[d] = 0;
							indices[d + 1]++;
						}
					}

					// is there another index after that?
					hasNext = false;
					for (int d=0; d<dimensions.length; d++) {
						if (indices[d] < dimensions[d] - 1) {
							hasNext = true;
							break;
						}
					}

					return indices;
				}
			};
		}
	}
}
