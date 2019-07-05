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
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.TupleTree;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.BigExp;

import java.math.MathContext;
import java.math.RoundingMode;


/** a matrix of Boltzmann-weighted energies */
public class ZMatrix extends TupleMatrixGeneric<BigExp> {

	// since we're outputting BigExp values, we only need about 16 decimal digits of precision
	private final BoltzmannCalculator bcalc = new BoltzmannCalculator(new MathContext(16, RoundingMode.HALF_UP));

	public ZMatrix(SimpleConfSpace confSpace) {
		super(confSpace);
	}

	private BigExp calc(double energy) {
		return new BigExp(bcalc.calcPrecise(energy));
	}

	public void set(EnergyMatrix emat) {

		// convert singles and pairs
		int n = getNumPos();
		for (int pos1=0; pos1<n; pos1++) {
			int m1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<m1; rc1++) {

				this.setOneBody(pos1, rc1, calc(emat.getOneBody(pos1, rc1)));

				for (int pos2=0; pos2<pos1; pos2++) {
					int m2 = getNumConfAtPos(pos2);
					for (int rc2=0; rc2<m2; rc2++) {

						this.setPairwise(pos1, rc1, pos2, rc2, calc(emat.getPairwise(pos1, rc1, pos2, rc2)));
					}
				}
			}
		}

		// convert higher-order tuples
		if (emat.hasHigherOrderTuples()) {

			// for each pair
			for (int pos1=0; pos1<n; pos1++) {
				int m1 = getNumConfAtPos(pos1);
				for (int rc1=0; rc1<m1; rc1++) {
					for (int pos2=0; pos2<pos1; pos2++) {
						int m2 = getNumConfAtPos(pos2);
						for (int rc2=0; rc2<m2; rc2++) {

							// convert all the higher-order tuples for that pair
							TupleTree<Double> tuples = emat.getHigherOrderTuples(pos1, rc1, pos2, rc2);
							if (tuples != null) {
								for (RCTuple tuple : tuples.makeTuplesList()) {
									this.setTuple(tuple, calc(tuples.get(tuple)));
								}
							}
						}
					}
				}
			}
		}
	}
}
