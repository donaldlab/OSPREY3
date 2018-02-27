package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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
		} else {
			return BigDecimal.valueOf(val);
		}
	}

	public static BigDecimal biggen(long val) {
		return new BigDecimal(val);
	}

	/**
	 * Tests for sameness of values,
	 * rather than BigDecimal.equals(), which tests sameness of representation
	 * (i.e., the same value can have multiple representations in BigDecimal)
	 */
	public static boolean isSameValue(BigDecimal a, BigDecimal b) {
		return a == b || a.compareTo(b) == 0;
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
		return !isInf(d) && d.compareTo(BigDecimal.ZERO) == 0;
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
		return !isInf(d) && d == BigNaN;
	}

	public static boolean isNaN(BigDecimal d) {
		return d == BigNaN;
	}

	/** return a < b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isLessThan(BigDecimal a, BigDecimal b) {
		if (a == BigNaN || b == BigNaN) {
			throw new IllegalArgumentException("can't compare NaN");
		}
		if (a == BigPositiveInfinity) {
			return false;
		} else if (a == BigNegativeInfinity) {
			if (b == BigNegativeInfinity) {
				return false;
			} else {
				return true;
			}
		} else {
			if (b == BigPositiveInfinity) {
				return true;
			} else if (b == BigNegativeInfinity) {
				return false;
			} else {
				return a.compareTo(b) < 0;
			}
		}
	}

	/** return a <= b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isLessThanOrEqual(BigDecimal a, BigDecimal b) {
		return isLessThan(a, b) || isSameValue(a, b);
	}

	/** return a > b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isGreaterThan(BigDecimal a, BigDecimal b) {
		if (a == BigNaN || b == BigNaN) {
			throw new IllegalArgumentException("can't compare NaN");
		}
		if (a == BigPositiveInfinity) {
			if (b == BigPositiveInfinity) {
				return false;
			} else {
				return true;
			}
		} else if (a == BigNegativeInfinity) {
			return false;
		} else {
			if (b == BigPositiveInfinity) {
				return false;
			} else if (b == BigNegativeInfinity) {
				return true;
			} else {
				return a.compareTo(b) > 0;
			}
		}
	}

	/** return a >= b, correctly handling -Inf, +Inf, and NaN */
	public static boolean isGreaterThanOrEqual(BigDecimal a, BigDecimal b) {
		return isGreaterThan(a, b) || isSameValue(a, b);
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
		if (a == BigNaN || b == BigNaN) {
			return BigNaN;
		}
		if (isZero(a) || isZero(b)) {
			return BigDecimal.ZERO;
		}
		if (a == BigPositiveInfinity) {
			if (b == BigNegativeInfinity) {
				return BigNaN;
			} else if (isPositive(b)) {
				return BigPositiveInfinity;
			} else {
				return BigNegativeInfinity;
			}
		} else if (a == BigNegativeInfinity) {
			if (b == BigPositiveInfinity) {
				return BigNaN;
			} else if (isPositive(b)) {
				return BigNegativeInfinity;
			} else {
				return BigPositiveInfinity;
			}
		} else if (isPositive(a)) {
			if (b == BigPositiveInfinity) {
				return BigPositiveInfinity;
			} else if (b == BigNegativeInfinity) {
				return BigNegativeInfinity;
			} else {
				return a.multiply(b, context);
			}
		} else {
			if (b == BigPositiveInfinity) {
				return BigNegativeInfinity;
			} else if (b == BigNegativeInfinity) {
				return BigPositiveInfinity;
			} else {
				return a.multiply(b, context);
			}
		}
	}
	// TODO: test me!!!

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
		return Math.log10(x.add(BigDecimal.ONE).doubleValue());
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
}
