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
