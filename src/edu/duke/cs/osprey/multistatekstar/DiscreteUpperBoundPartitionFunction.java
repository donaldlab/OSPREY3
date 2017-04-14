package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator.Async;

import java.math.BigDecimal;
import java.math.RoundingMode;

import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 * Computes a 1+epsilon approximation to the partition function
 */
public class DiscreteUpperBoundPartitionFunction extends DiscretePartitionFunction {

	public DiscreteUpperBoundPartitionFunction(EnergyMatrix emat, PruningMatrix pmat,
			ConfSearchFactory confSearchFactory, Async ecalc) {
		super(emat, pmat, confSearchFactory, ecalc);
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
		if(status == Status.Estimated) {
			values.qstar = values.qstar.multiply(new BigDecimal(1.0+getEffectiveEpsilon()));
		}
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