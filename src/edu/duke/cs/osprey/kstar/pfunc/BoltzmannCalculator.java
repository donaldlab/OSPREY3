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
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.MathTools.BigDecimalBounds;
import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import static edu.duke.cs.osprey.tools.Log.log;


public class BoltzmannCalculator {

	public static double constRT = PoissonBoltzmannEnergy.constRT;

	public final MathContext mathContext;
	public final ExpFunction e;

	public BoltzmannCalculator(MathContext mathContext) {
		this.mathContext = mathContext;
		this.e = new ExpFunction(mathContext);
	}
	
	public BigDecimal calc(double energy) {
		return e.exp(-energy/constRT);
	}

	public double freeEnergy(BigDecimal z) {
		return -constRT*e.log(z).doubleValue();
	}

	public BigDecimal calcPrecise(double e) {
		return exp(-e/constRT);
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
		return -constRT*ln(z);
	}

	public double freeEnergyPrecise(BigExp z) {
		return freeEnergyPrecise(z.toBigDecimal(mathContext));
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
