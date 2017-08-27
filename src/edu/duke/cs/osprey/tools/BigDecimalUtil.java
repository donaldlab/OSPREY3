package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * BigDecimal utilities.
 *
 * @author Valentyn Kolesnikov
 * @version $Revision$ $Date$
 */
public final class BigDecimalUtil {
	private static final int SCALE = 64;
	public static long ITER = 1000;
	public static MathContext context = new MathContext( 100 );
	private static final int ROUNDING_MODE = BigDecimal.ROUND_HALF_EVEN;
	public static BigDecimal PI_DIV_180
	= new BigDecimal("3.1415926535897932384626433832795").divide(BigDecimal.valueOf(180), 32, BigDecimal.ROUND_HALF_UP);
	public static BigDecimal PI_DIV_200
	= new BigDecimal("3.1415926535897932384626433832795").divide(BigDecimal.valueOf(200), 32, BigDecimal.ROUND_HALF_UP);
	public static BigDecimal EPS = BigDecimal.ONE.scaleByPowerOfTen(-100);

	private BigDecimalUtil() {
	}

	/**
	 * Compute the square root of x to a given scale, x &gt;= 0.
	 * Use Newton's algorithm.
	 * @param x the value of x
	 * @return the result value
	 */
	public static BigDecimal sqrt(BigDecimal x) {
		// Check that x >= 0.
		if (x.signum() < 0) {
			throw new ArithmeticException("x < 0");
		}

		// n = x*(10^(2*SCALE))
		BigInteger n = x.movePointRight(SCALE << 1).toBigInteger();

		// The first approximation is the upper half of n.
		int bits = (n.bitLength() + 1) >> 1;
		BigInteger ix = n.shiftRight(bits);
		BigInteger ixPrev;

		// Loop until the approximations converge
		// (two successive approximations are equal after rounding).
		do {
			ixPrev = ix;

			// x = (x + n/x)/2
			ix = ix.add(n.divide(ix)).shiftRight(1);

			Thread.yield();
		} while (ix.compareTo(ixPrev) != 0);

		return new BigDecimal(ix, SCALE);
	}

	/**
	 * Compute the integral root of x to a given scale, x &gt;= 0.
	 * Use Newton's algorithm.
	 * @param x the value of x
	 * @param index the integral root value
	 * @param scale the desired scale of the result
	 * @return the result value
	 */
	public static BigDecimal intRoot(BigDecimal x, long index,
			int scale)
	{
		// Check that x >= 0.
		if (x.signum() < 0) {
			throw new IllegalArgumentException("x < 0");
		}

		int        sp1 = scale + 1;
		BigDecimal n   = x;
		BigDecimal i   = BigDecimal.valueOf(index);
		BigDecimal im1 = BigDecimal.valueOf(index-1);
		BigDecimal tolerance = BigDecimal.valueOf(5)
				.movePointLeft(sp1);
		BigDecimal xPrev;

		// The initial approximation is x/index.
		x = x.divide(i, scale, BigDecimal.ROUND_HALF_EVEN);

		// Loop until the approximations converge
		// (two successive approximations are equal after rounding).
		do {
			// x^(index-1)
			BigDecimal xToIm1 = intPower(x, index-1, sp1);

			// x^index
			BigDecimal xToI =
					x.multiply(xToIm1)
					.setScale(sp1, BigDecimal.ROUND_HALF_EVEN);

			// n + (index-1)*(x^index)
			BigDecimal numerator =
					n.add(im1.multiply(xToI))
					.setScale(sp1, BigDecimal.ROUND_HALF_EVEN);

			// (index*(x^(index-1))
			BigDecimal denominator =
					i.multiply(xToIm1)
					.setScale(sp1, BigDecimal.ROUND_HALF_EVEN);

			// x = (n + (index-1)*(x^index)) / (index*(x^(index-1)))
			xPrev = x;
			x = numerator
					.divide(denominator, sp1, BigDecimal.ROUND_DOWN);

			Thread.yield();
		} while (x.subtract(xPrev).abs().compareTo(tolerance) > 0);

		return x;
	}

	/**
	 * Compute the natural logarithm of x to a given scale, x &gt; 0.
	 * @param x the value
	 * @param scale the scale
	 * @return the result
	 */
	public static BigDecimal ln(BigDecimal x, int scale) {
		// Check that x > 0.
		if (x.signum() <= 0) {
			throw new IllegalArgumentException("x <= 0");
		}

		// The number of digits to the left of the decimal point.
		int magnitude = x.toString().length() - x.scale() - 1;

		if (magnitude < 3) {
			return lnNewton(x, scale);
		}

		// Compute magnitude*ln(x^(1/magnitude)).
		else {

			// x^(1/magnitude)
			BigDecimal root = intRoot(x, magnitude, scale);

			// ln(x^(1/magnitude))
			BigDecimal lnRoot = lnNewton(root, scale);

			// magnitude*ln(x^(1/magnitude))
			return BigDecimal.valueOf(magnitude).multiply(lnRoot)
					.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
		}
	}

	/**
	 * Compute the natural logarithm of x to a given scale, x > 0.
	 * Use Newton's algorithm.
	 */
	private static BigDecimal lnNewton(BigDecimal x, int scale)
	{
		int        sp1 = scale + 1;
		BigDecimal n   = x;
		BigDecimal term;

		// Convergence tolerance = 5*(10^-(scale+1))
		BigDecimal tolerance = BigDecimal.valueOf(5)
				.movePointLeft(sp1);

		// Loop until the approximations converge
		// (two successive approximations are within the tolerance).
		do {

			// e^x
			BigDecimal eToX = exp(x, sp1);

			// (e^x - n)/e^x
			term = eToX.subtract(n)
					.divide(eToX, sp1, BigDecimal.ROUND_DOWN);

			// x - (e^x - n)/e^x
			x = x.subtract(term);

			Thread.yield();
		} while (term.compareTo(tolerance) > 0);

		return x.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
	}

	public static BigDecimal cosine(BigDecimal x) {

		BigDecimal currentValue = BigDecimal.ONE;
		BigDecimal lastVal      = currentValue.add(BigDecimal.ONE);
		BigDecimal xSquared     = x.multiply(x);
		BigDecimal numerator    = BigDecimal.ONE;
		BigDecimal denominator  = BigDecimal.ONE;
		int        i            = 0;

		while (lastVal.compareTo(currentValue) != 0) {
			lastVal = currentValue;

			int z = 2 * i + 2;

			denominator = denominator.multiply(BigDecimal.valueOf(z));
			denominator = denominator.multiply(BigDecimal.valueOf(z - 1));
			numerator   = numerator.multiply(xSquared);

			BigDecimal term = numerator.divide(denominator, SCALE + 5, ROUNDING_MODE);

			if (i % 2 == 0) {
				currentValue = currentValue.subtract(term);
			} else {
				currentValue = currentValue.add(term);
			}
			i++;
		}

		return currentValue;
	}

	public static BigDecimal sine(BigDecimal x) {
		BigDecimal lastVal      = x.add(BigDecimal.ONE);
		BigDecimal currentValue = x;
		BigDecimal xSquared     = x.multiply(x);
		BigDecimal numerator    = x;
		BigDecimal denominator  = BigDecimal.ONE;
		int        i            = 0;

		while (lastVal.compareTo(currentValue) != 0) {
			lastVal = currentValue;

			int z = 2 * i + 3;

			denominator = denominator.multiply(BigDecimal.valueOf(z));
			denominator = denominator.multiply(BigDecimal.valueOf(z - 1));
			numerator   = numerator.multiply(xSquared);

			BigDecimal term = numerator.divide(denominator, SCALE + 5, ROUNDING_MODE);

			if (i % 2 == 0) {
				currentValue = currentValue.subtract(term);
			} else {
				currentValue = currentValue.add(term);
			}

			i++;
		}
		return currentValue;
	}

	public static BigDecimal tangent(BigDecimal x) {

		BigDecimal sin = sine(x);
		BigDecimal cos = cosine(x);

		return sin.divide(cos, SCALE, BigDecimal.ROUND_HALF_UP);
	}

	public static BigDecimal log10(BigDecimal b) {
		final int NUM_OF_DIGITS = SCALE + 2;
		// need to add one to get the right number of dp
		//  and then add one again to get the next number
		//  so I can round it correctly.

		MathContext mc = new MathContext(NUM_OF_DIGITS, RoundingMode.HALF_EVEN);
		//special conditions:
		// log(-x) -> exception
		// log(1) == 0 exactly;
		// log of a number lessthan one = -log(1/x)
		if(b.signum() <= 0) {
			throw new ArithmeticException("log of a negative number! (or zero)");
		}
		else if(b.compareTo(BigDecimal.ONE) == 0) {
			return BigDecimal.ZERO;
		} else if(b.compareTo(BigDecimal.ONE) < 0) {
			return (log10((BigDecimal.ONE).divide(b,mc))).negate();
		}

		StringBuilder sb = new StringBuilder();
		//number of digits on the left of the decimal point
		int leftDigits = b.precision() - b.scale();

		//so, the first digits of the log10 are:
		sb.append(leftDigits - 1).append(".");

		//this is the algorithm outlined in the webpage
		int n = 0;
		while (n < NUM_OF_DIGITS) {
			b = (b.movePointLeft(leftDigits - 1)).pow(10, mc);
			leftDigits = b.precision() - b.scale();
			sb.append(leftDigits - 1);
			n++;
		}

		BigDecimal ans = new BigDecimal(sb.toString());

		//Round the number to the correct number of decimal places.
		ans = ans.round(new MathContext(ans.precision() - ans.scale() + SCALE, RoundingMode.HALF_EVEN));
		return ans;
	}

	public static BigDecimal cuberoot(BigDecimal b) {
		// Specify a math context with 40 digits of precision.

		MathContext mc = new MathContext(40);

		BigDecimal x = new BigDecimal("1", mc);

		// Search for the cube root via the Newton-Raphson loop. Output each // successive iteration's value.

		for (int i = 0; i < ITER; i++) {
			x = x.subtract(x.pow(3, mc).subtract(b, mc).divide(new BigDecimal("3", mc).multiply(x.pow(2, mc), mc), mc), mc);
		}
		return x;
	}

	public static BigDecimal pow(BigDecimal savedValue, BigDecimal value) {
		BigDecimal result = null;
		result = exp(ln(savedValue, 32).multiply(value), 32);
		return result;
	}

	/**
	 * Compute x^exponent to a given scale.  Uses the same
	 * algorithm as class numbercruncher.mathutils.IntPower.
	 * @param x the value x
	 * @param exponent the exponent value
	 * @param scale the desired scale of the result
	 * @return the result value
	 */
	public static BigDecimal intPower(BigDecimal x, long exponent,
			int scale) {
		// If the exponent is negative, compute 1/(x^-exponent).
		if (exponent < 0) {
			return BigDecimal.valueOf(1)
					.divide(intPower(x, -exponent, scale), scale,
							BigDecimal.ROUND_HALF_EVEN);
		}

		BigDecimal power = BigDecimal.valueOf(1);

		// Loop to compute value^exponent.
		while (exponent > 0) {

			// Is the rightmost bit a 1?
			if ((exponent & 1) == 1) {
				power = power.multiply(x)
						.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
			}

			// Square x and shift exponent 1 bit to the right.
			x = x.multiply(x)
					.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
			exponent >>= 1;

		Thread.yield();
		}

		return power;
	}

	/**
	 * Compute e^x to a given scale.
	 * Break x into its whole and fraction parts and
	 * compute (e^(1 + fraction/whole))^whole using Taylor's formula.
	 * @param x the value of x
	 * @param scale the desired scale of the result
	 * @return the result value
	 */
	public static BigDecimal exp(BigDecimal x, int scale) {
		// e^0 = 1
		if (x.signum() == 0) {
			return BigDecimal.valueOf(1);
		}

		// If x is negative, return 1/(e^-x).
		else if (x.signum() == -1) {
			return BigDecimal.valueOf(1)
					.divide(exp(x.negate(), scale), scale,
							BigDecimal.ROUND_HALF_EVEN);
		}

		// Compute the whole part of x.
		BigDecimal xWhole = x.setScale(0, BigDecimal.ROUND_DOWN);

		// If there isn't a whole part, compute and return e^x.
		if (xWhole.signum() == 0) {
			return expTaylor(x, scale);
		}

		// Compute the fraction part of x.
		BigDecimal xFraction = x.subtract(xWhole);

		// z = 1 + fraction/whole
		BigDecimal z = BigDecimal.valueOf(1)
				.add(xFraction.divide(
						xWhole, scale,
						BigDecimal.ROUND_HALF_EVEN));

		// t = e^z
		BigDecimal t = expTaylor(z, scale);

		BigDecimal maxLong = BigDecimal.valueOf(Long.MAX_VALUE);
		BigDecimal result  = BigDecimal.valueOf(1);

		// Compute and return t^whole using intPower().
		// If whole > Long.MAX_VALUE, then first compute products
		// of e^Long.MAX_VALUE.
		while (xWhole.compareTo(maxLong) >= 0) {
			result = result.multiply(
					intPower(t, Long.MAX_VALUE, scale))
					.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
			xWhole = xWhole.subtract(maxLong);

			Thread.yield();
		}
		return result.multiply(intPower(t, xWhole.longValue(), scale))
				.setScale(scale, BigDecimal.ROUND_HALF_EVEN);
	}

	/**
	 * Compute e^x to a given scale by the Taylor series.
	 * @param x the value of x
	 * @param scale the desired scale of the result
	 * @return the result value
	 */
	private static BigDecimal expTaylor(BigDecimal x, int scale)
	{
		BigDecimal factorial = BigDecimal.valueOf(1);
		BigDecimal xPower    = x;
		BigDecimal sumPrev;

		// 1 + x
		BigDecimal sum  = x.add(BigDecimal.valueOf(1));

		// Loop until the sums converge
		// (two successive sums are equal after rounding).
		int i = 2;
		do {
			// x^i
			xPower = xPower.multiply(x)
					.setScale(scale, BigDecimal.ROUND_HALF_EVEN);

			// i!
			factorial = factorial.multiply(BigDecimal.valueOf(i));

			// x^i/i!
			BigDecimal term = xPower
					.divide(factorial, scale,
							BigDecimal.ROUND_HALF_EVEN);

			// sum = sum + x^i/i!
			sumPrev = sum;
			sum = sum.add(term);

			++i;
			Thread.yield();
		} while (sum.compareTo(sumPrev) != 0);

		return sum;
	}

	public static BigDecimal asin(BigDecimal val) {
		return BigDecimal.valueOf(Math.asin(val.doubleValue()));
	}

	public static BigDecimal acos(BigDecimal val) {
		return BigDecimal.valueOf(Math.acos(val.doubleValue()));
	}

	public static BigDecimal atan(BigDecimal val) {
		return BigDecimal.valueOf(Math.atan(val.doubleValue()));
	}
}
