package edu.duke.cs.osprey.tools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;


/**
 * A ridiculously simplistic and naive way to have floating point values with larger exponents than +-300 ish
 * but also have fast performance for multiply and max
 * (BigDecimal is relatively slow)
 *
 * The BigExp number type consists of a floating point component (fp), and an exponent component (exp).
 * fp is conditioned so its decimal exponent lies within [-100,100], with any overflow saved in exp.
 *
 * Benchmarking indicates BigExp is roughly 30x faster than BigDecimal multiplication
 * and roughly 2.4x slower than native double multiplication
 */
public class BigExp implements Comparable<BigExp>, IOable {

	// since we're using doubles for the significand, we only ever need ~16 digits of precision
	public static final MathContext mathContext = new MathContext(16, RoundingMode.HALF_UP);

	public double fp;
	public int exp;

	public BigExp(double fp, int exp) {
		this.fp = fp;
		this.exp = exp;
	}

	public BigExp(double fp) {
		set(fp);
	}

	public BigExp(BigExp other) {
		set(other);
	}

	public BigExp(BigDecimal bd) {
		set(bd);
	}

	public BigExp(BigInteger bi) {
		set(bi);
	}

	public void set(double fp) {
		this.fp = fp;
		this.exp = 0;
		normalize(false);
	}

	public void set(double fp, int exp) {
		this.fp = fp;
		this.exp = exp;
	}

	public void set(BigExp other) {
		set(other.fp, other.exp);
	}

	public void set(BigDecimal bd) {

		// not super efficient, but we don't need this part to be crazy fast
		if (MathTools.isFinite(bd)) {
			String[] str = String.format("%.16e", bd).split("e");
			fp = Double.parseDouble(str[0]);
			exp = Integer.parseInt(str[1]);
			// no normalization needed
			//assert ((fp < 1e1 && fp >= 1e-1) || (fp >= -1e-1 && fp > -1e1));
		} else {
			fp = bd.doubleValue();
			exp = 0;
		}
	}

	public void set(BigInteger bi) {
		set(new BigDecimal(bi));
	}

	public void normalize(boolean fully) {

		if (!Double.isFinite(fp) || fp == 0.0) {
			exp = 0;
			return;
		}

		// this is super naive, but it works and it's actually really fast! =D

		if (fp >= 1) {

			while (fp > 1e100) {
				fp *= 1e-100;
				exp += 100;
			}

			if (fully) {

				while (fp > 1e10) {
					fp *= 1e-10;
					exp += 10;
				}

				while (fp > 1e1) {
					fp *= 1e-1;
					exp += 1;
				}
			}

		} else if (fp > 0) {

			while (fp < 1e-100) {
				fp *= 1e100;
				exp -= 100;
			}

			if (fully) {

				while (fp < 1e-10) {
					fp *= 1e10;
					exp -= 10;
				}

				while (fp < 1e-1) {
					fp *= 1e1;
					exp -= 1;
				}
			}

		} else if (fp > -1) {

			while (fp > -1e-100) {
				fp *= 1e100;
				exp -= 100;
			}

			if (fully) {

				while (fp > -1e-10) {
					fp *= 1e10;
					exp -= 10;
				}

				while (fp > -1e-1) {
					fp *= 1e1;
					exp -= 1;
				}
			}

		} else { // if (fp <= -1)

			while (fp < -1e100) {
				fp *= 1e-100;
				exp += 100;
			}

			if (fully) {

				while (fp < -1e10) {
					fp *= 1e-10;
					exp += 10;
				}

				while (fp < -1e1) {
					fp *= 1e-1;
					exp += 1;
				}
			}
		}
	}

	public void mult(BigExp other) {
		this.fp *= other.fp;
		this.exp += other.exp;
		normalize(false);
	}

	public void mult(double other) {
		this.fp *= other;
		normalize(false);
	}

	public void mult(int other) {
		this.fp *= other;
		normalize(false);
	}

	public void mult(long other) {
		this.fp *= other;
		normalize(false);
	}

	public void div(BigExp other) {
		this.fp /= other.fp;
		this.exp -= other.exp;
		normalize(false);
	}

	public void div(double other) {
		this.fp /= other;
		normalize(false);
	}

	public void div(int other) {
		this.fp /= other;
		normalize(false);
	}

	public void div(long other) {
		this.fp /= other;
		normalize(false);
	}

	public void add(BigExp other) {
		set(this.toBigDecimal().add(other.toBigDecimal(), mathContext));
	}

	public void sub(BigExp other) {
		set(this.toBigDecimal().subtract(other.toBigDecimal(), mathContext));
	}

	public void abs() {
		fp = Math.abs(fp);
	}

	public void min(BigExp other) {
		if (compareTo(other) > 0) {
			set(other);
		}
	}

	public void max(BigExp other) {
		if (compareTo(other) < 0) {
			set(other);
		}
	}

	@Override
	public int compareTo(BigExp other) {

		// since we're comparing exponents directly, we need to use full normalization
		this.normalize(true);
		other.normalize(true);

		double thisSign = this.sign();
		double otherSign = other.sign();
		if (thisSign > otherSign) {
			return 1;
		} else if (thisSign < otherSign) {
			return -1;
		} else {

			int sign = (int)thisSign;
			if (this.exp > other.exp) {
				return sign;
			} else if (this.exp < other.exp) {
				return -sign;
			}

			return Double.compare(this.fp, other.fp);
		}
	}

	public boolean greaterThan(BigExp other) {
		return compareTo(other) > 0;
	}

	public boolean greaterThanOrEqual(BigExp other) {
		return compareTo(other) >= 0;
	}

	public boolean lessThan(BigExp other) {
		return compareTo(other) < 0;
	}

	public boolean lessThanOrEqual(BigExp other) {
		return compareTo(other) <= 0;
	}

	public boolean isNaN() {
		return Double.isNaN(fp);
	}

	public boolean isFinite() {
		return Double.isFinite(fp);
	}

	public double sign() {
		return Math.signum(fp);
	}

	public boolean isPositive() {
		return fp > 0;
	}

	public boolean isNegative() {
		return fp < 0;
	}

	public BigDecimal toBigDecimal() {
		return toBigDecimal(MathContext.UNLIMITED);
	}

	public BigDecimal toBigDecimal(MathContext mc) {
		if (Double.isNaN(fp)) {
			return MathTools.BigNaN;
		} else if (fp == Double.POSITIVE_INFINITY) {
			return MathTools.BigPositiveInfinity;
		} else if (fp == Double.NEGATIVE_INFINITY) {
			return MathTools.BigNegativeInfinity;
		} else {

			// fully normalize before converting to big decimal,
			// to avoid some roundoff error when different normalization states get convert to big decimals
			normalize(true);
			BigDecimal bd = BigDecimal.valueOf(fp);
			return new BigDecimal(bd.unscaledValue(), bd.scale() - exp, mc);
		}
	}

	public double toDouble() {
		// TODO: optimize this?
		return toBigDecimal(mathContext).doubleValue();
	}

	@Override
	public int hashCode() {
		return HashCalculator.combineHashes(
			Double.hashCode(fp),
			Integer.hashCode(exp)
		);
	}

	@Override
	public boolean equals(Object other) {
		return other instanceof BigExp && equals((BigExp)other);
	}

	public boolean equals(BigExp other) {
		return this.compareTo(other) == 0;
	}

	@Override
	public String toString() {
		return toString(6);
	}

	public String toString(int precision) {
		String[] parts = String.format("%." + precision + "e", fp).split("e");
		if (parts.length == 2) {
			return parts[0] + "e" + (Integer.parseInt(parts[1]) + exp);
		} else {
			return parts[0];
		}
	}

	public static final int NumBytes = 12;

	@Override
	public void writeTo(DataOutput out)
	throws IOException {
		out.writeDouble(fp);
		out.writeInt(exp);
	}

	@Override
	public void readFrom(DataInput in)
	throws IOException {
		fp = in.readDouble();
		exp = in.readInt();
	}

	public static BigExp read(DataInput in)
	throws IOException {
		BigExp out = new BigExp(0.0, 0);
		out.readFrom(in);
		return out;
	}


	public static class Bounds {

		public BigExp lower;
		public BigExp upper;

		public Bounds(BigExp lower, BigExp upper) {
			this.lower = lower;
			this.upper = upper;
		}

		public Bounds(BigExp val) {
			this(val, val);
		}

		public Bounds() {
			this(new BigExp(Double.NEGATIVE_INFINITY), new BigExp(Double.POSITIVE_INFINITY));
		}

		public boolean isValid() {
			return lower.lessThanOrEqual(upper);
		}

		public boolean contains(BigExp query) {
			return query.greaterThanOrEqual(lower)
				&& query.lessThanOrEqual(upper);
		}

		public boolean contains(BigExp.Bounds query) {
			return query.lower.greaterThanOrEqual(lower)
				&& query.upper.lessThanOrEqual(upper);
		}

		@Override
		public String toString() {
			return String.format("[%s,%s]", lower, upper);
		}
	}
}
