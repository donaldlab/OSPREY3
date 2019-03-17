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

package edu.duke.cs.osprey.sofea;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.lute.LUTEConfEnergyCalculator;

import java.math.BigDecimal;
import java.math.MathContext;


public class BoltzmannLute {

	public final LUTEConfEnergyCalculator luteEcalc;

	public final BigDecimal factor;
	private final BigDecimal[] values;

	public BoltzmannLute(LUTEConfEnergyCalculator luteEcalc, MathContext mathContext) {

		this.luteEcalc = luteEcalc;
		this.values = new BigDecimal[luteEcalc.tuples.size()];

		// pre-calculate all the boltzmann-weighted values
		BoltzmannCalculator bcalc = new BoltzmannCalculator(mathContext);
		this.factor = bcalc.calcPrecise(luteEcalc.state.tupleEnergyOffset);
		for (int t=0; t<luteEcalc.tuples.size(); t++) {
			this.values[t] = bcalc.calcPrecise(luteEcalc.state.tupleEnergies[t]);
		}
	}

	public boolean has(int pos, int rc) {
		return luteEcalc.hasTuple(pos, rc);
	}

	public boolean has(int pos1, int rc1, int pos2, int rc2) {
		return luteEcalc.hasTuple(pos1, rc1, pos2, rc2);
	}

	public BigDecimal get(int pos, int rc) {
		return get(luteEcalc.tuples.getIndex(pos, rc));
	}

	public BigDecimal get(int pos1, int rc1, int pos2, int rc2) {
		return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2));
	}

	public BigDecimal get(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return get(luteEcalc.tuples.getIndex(pos1, rc1, pos2, rc2, pos3, rc3));
	}

	public BigDecimal get(RCTuple tuple) {
		return get(luteEcalc.tuples.getIndex(tuple));
	}

	private BigDecimal get(Integer t) {
		if (t == null) {
			return null;
		} else {
			return values[t];
		}
	}
}
