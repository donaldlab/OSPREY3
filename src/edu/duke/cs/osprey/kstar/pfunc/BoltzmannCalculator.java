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
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;

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
			return BigDecimalMath.exp(new BigDecimal(e), mathContext);
		}
	}

	public double freeEnergyPrecise(BigDecimal z) {
		return -constRT*ln(z);
	}

	public double ln(BigDecimal z) {
		if (MathTools.isNaN(z) || MathTools.isNegative(z)) {
			return Double.NaN;
		} else if (MathTools.isZero(z)) {
			return Double.NEGATIVE_INFINITY;
		} else if (z == MathTools.BigPositiveInfinity) {
			return Double.POSITIVE_INFINITY;
		} else {
			return BigDecimalMath.log(z, mathContext).doubleValue();
		}
	}

	public double ln1p(BigDecimal z) {
		return ln(z.add(BigDecimal.ONE));
	}
}
