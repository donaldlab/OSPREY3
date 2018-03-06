package edu.duke.cs.osprey.tools;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;

public class BigMath {

	public final MathContext context;

	private BigDecimal d;

	public BigMath(MathContext context) {
		this.context = context;
	}

	public BigMath(int precision) {
		this(new MathContext(precision, RoundingMode.HALF_UP));
	}

	public BigDecimal get() {
		return d;
	}
	public BigMath set(BigDecimal val) {
		d = val;
		return this;
	}
	public BigMath set(BigInteger val) {
		return set(new BigDecimal(val));
	}
	public BigMath set(double val) {
		return set(MathTools.biggen(val));
	}
	public BigMath set(long val) {
		return set(MathTools.biggen(val));
	}

	public BigMath add(BigDecimal other) {
		d = MathTools.bigAdd(d, other, context);
		return this;
	}
	public BigMath add(BigInteger other) {
		return add(new BigDecimal(other));
	}
	public BigMath add(double val) {
		return add(MathTools.biggen(val));
	}
	public BigMath add(long val) {
		return add(MathTools.biggen(val));
	}

	public BigMath sub(BigDecimal other) {
		d = MathTools.bigSubtract(d, other, context);
		return this;
	}
	public BigMath sub(BigInteger other) {
		return sub(new BigDecimal(other));
	}
	public BigMath sub(double val) {
		return sub(MathTools.biggen(val));
	}
	public BigMath sub(long val) {
		return sub(MathTools.biggen(val));
	}

	public BigMath mult(BigDecimal other) {
		d = MathTools.bigMultiply(d, other,context);
		return this;
	}
	public BigMath mult(BigInteger other) {
		return mult(new BigDecimal(other));
	}
	public BigMath mult(double val) {
		return mult(MathTools.biggen(val));
	}
	public BigMath mult(long val) {
		return mult(MathTools.biggen(val));
	}

	public BigMath div(BigDecimal other) {
		d = MathTools.bigDivide(d, other, context);
		return this;
	}
	public BigMath div(BigInteger other) {
		return div(new BigDecimal(other));
	}
	public BigMath div(double val) {
		return div(MathTools.biggen(val));
	}
	public BigMath div(long val) {
		return div(MathTools.biggen(val));
	}
}
