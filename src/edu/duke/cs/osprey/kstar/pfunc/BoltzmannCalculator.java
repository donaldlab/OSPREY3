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

package edu.duke.cs.osprey.kstar.pfunc;

import java.math.BigDecimal;
import java.math.MathContext;

import ch.obermuhlner.math.big.BigDecimalMath;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import static edu.duke.cs.osprey.tools.Log.log;


public class BoltzmannCalculator {

	public static final double RClassic = 1.9891/1000.0;
	public static final double RPrecise = 1.98720425864083e-3; // from wikipedia

	public static final double TAlmostFrozen = 275; // about 35.33 F
	public static final double TRoom = 298.15; // exactly 25 C, 77 F
	public static final double TBody = 310; // about 98.33 F
	public static final double TClassic = TRoom;

	public static class Conditions {

		/** gas constant, in kcal/K/mol */
		public final double R;

		/** temperature in Kelvins */
		public final double T;

		public Conditions(double R, double T) {
			this.R = R;
			this.T = T;
		}

		public Conditions(double T) {
			this(RPrecise, T);
		}

		/** The time-tested consitions used in Osprey versions of years past */
		public static final Conditions Classic = new Conditions(RClassic, TClassic);

		public static final Conditions AlmostFrozen = new Conditions(TAlmostFrozen);
		public static final Conditions Room = new Conditions(TRoom);
		public static final Conditions Body = new Conditions(TBody);
	}

	public final MathContext mathContext;
	public final ExpFunction e;

	public final double R;

	public final double T;

	private final double RT;

	public BoltzmannCalculator(MathContext mathContext) {
		this(mathContext, RClassic, TClassic);
	}

	public BoltzmannCalculator(MathContext mathContext, double R, double T) {
		this.mathContext = mathContext;
		this.e = new ExpFunction(mathContext);
		this.R = R;
		this.T = T;
		this.RT = R*T;
	}

	public BoltzmannCalculator(MathContext mathContext, Conditions conditions) {
		this(mathContext, conditions.R, conditions.T);
	}
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/RT);
	}

	public double freeEnergy(BigDecimal z) {
		return -RT*e.log(z).doubleValue();
	}

	public BigDecimal calcPrecise(double e) {
		return exp(-e/RT);
	}

	public BigDecimal exp(double e) {
		if (Double.isNaN(e)) {
			return MathTools.BigNaN;
		} else if (e == Double.NEGATIVE_INFINITY) {
			return BigDecimal.ZERO;
		} else if (e == Double.POSITIVE_INFINITY) {
			return MathTools.BigPositiveInfinity;
		} else {

			// if the energy is too large, just return 0
			// NOTE: this logic checks the exact conditions in BigDecimalMath.exp() that triggers an overflow exception
			BigDecimal d = new BigDecimal(e);
			BigDecimal i = BigDecimalMath.integralPart(d);
			if (i.signum() != 0) {
				long num = i.longValueExact();
				if ((int)num != num) {
					if (e > 0) {
						//log("WARN: assuming exp of huge positive value (%.2e) is just +inf", e);
						return MathTools.BigPositiveInfinity;
					} else {
						//log("WARN: assuming exp of huge negative value (%.2e) is just 0", e);
						return BigDecimal.ZERO;
					}
				}
			}

			return BigDecimalMath.exp(d, mathContext);
		}
	}

	public double freeEnergyPrecise(BigDecimal z) {
		return -RT*ln(z);
	}

	public double freeEnergyPrecise(BigExp z) {
		return -RT*z.ln();
	}

	public DoubleBounds freeEnergyPrecise(BigDecimalBounds z) {
		// remember to swap the bounds, since computing free energy negates the value
		return new DoubleBounds(
			freeEnergyPrecise(z.upper),
			freeEnergyPrecise(z.lower)
		);
	}

	public DoubleBounds freeEnergyPrecise(BigExp.Bounds z) {
		// remember to swap the bounds, since computing free energy negates the value
		return new DoubleBounds(
			freeEnergyPrecise(z.upper),
			freeEnergyPrecise(z.lower)
		);
	}

	public double ln(BigDecimal z) {
		if (MathTools.isNaN(z) || MathTools.isNegative(z)) {
			return Double.NaN;
		} else if (MathTools.isZero(z)) {
			return Double.NEGATIVE_INFINITY;
		} else if (z == MathTools.BigPositiveInfinity) {
			return Double.POSITIVE_INFINITY;
		} else {

			// HACKHACK: values with really huge exponents take forever to calculate
			// if the exponent is too big, just cheat and say the answer is inf
			if (z.scale() > 10000) {
				log("WARN: assuming ln of tiny value (%.2e) is just -inf", z);
				return Double.NEGATIVE_INFINITY;
			} else if (z.scale() < -10000) {
				log("WARN: assuming ln of huge value (%.2e) is just +inf", z);
				return Double.POSITIVE_INFINITY;
			}

			return BigDecimalMath.log(z, mathContext).doubleValue();
		}
	}

	public double ln(BigExp z) {
		return ln(z.toBigDecimal(mathContext));
	}

	/**
	 * computes the following:
	 *  for z == 0: 0
	 *  for z > 0: ln(1+z)
	 *  for z < 0: -ln(1-z)
	 */
	public double ln1p(BigDecimal z) {
		double d;
		if (MathTools.isInf(z) || MathTools.isNaN(z) || MathTools.isZero(z)) {
			d = z.doubleValue();
		} else if (MathTools.isPositive(z)) {
			d = ln(z.add(BigDecimal.ONE, mathContext));
		} else {
			d = -ln(MathTools.bigNegate(z).add(BigDecimal.ONE, mathContext));
		}
		return d;
	}

	public double ln1p(BigExp z) {
		return ln1p(z.toBigDecimal(mathContext));
	}
}
