package edu.duke.cs.osprey.sofea;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

import java.math.BigDecimal;

/** a matrix of Boltzmann-weighted energies */
public class ZMatrix extends TupleMatrixGeneric<BigDecimal> {

	public ZMatrix(SimpleConfSpace confSpace) {
		super(confSpace);
	}

	public void set(EnergyMatrix emat, BoltzmannCalculator bcalc) {
		int n = getNumPos();
		for (int pos1=0; pos1<n; pos1++) {
			int m1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<m1; rc1++) {

				this.setOneBody(pos1, rc1, bcalc.calcPrecise(emat.getOneBody(pos1, rc1)));

				for (int pos2=0; pos2<pos1; pos2++) {
					int m2 = getNumConfAtPos(pos2);
					for (int rc2=0; rc2<m2; rc2++) {

						this.setPairwise(pos1, rc1, pos2, rc2, bcalc.calcPrecise(emat.getPairwise(pos1, rc1, pos2, rc2)));
					}
				}
			}
		}
	}
}
