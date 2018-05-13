/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator.Async;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Computes a 1+epsilon approximation to the partition function
 */
public class PartitionFunctionDiscreteUppperBound extends PartitionFunctionDiscrete {

	public PartitionFunctionDiscreteUppperBound(
			EnergyMatrix emat, 
			PruningMatrix pmat,
			PruningMatrix invmat,
			ConfSearchFactory confSearchFactory, 
			Async ecalc
			) {
		super(emat, pmat, invmat, confSearchFactory, ecalc);
	}

	@Override
	public void compute(int maxNumConfs) {
		super.compute(maxNumConfs);
		if(status == Status.Estimated) {
			values.qstar = values.qstar.multiply(new BigDecimal(1.0+getEffectiveEpsilon()));
		}
	}

	@Override
	public void compute(BigDecimal targetScoreWeights) {
		super.compute(targetScoreWeights);
		values.qstar = values.qstar.multiply(new BigDecimal(1.0+getEffectiveEpsilon()));
	}

	@Override
	protected double getEffectiveEpsilon() {
		BigDecimal s = values.qprime.add(values.pstar);
		BigDecimal q = s.add(values.qstar);

		if (q.compareTo(BigDecimal.ZERO) == 0) {
			// this is really bad... it should never happen
			// it probably means there are zero conformations in the tree
			return Double.NaN;
		}

		return s.divide(values.qstar, RoundingMode.HALF_UP).doubleValue();
	}

}
