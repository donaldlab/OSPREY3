package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public class KStarScoreUpperBound extends KStarScoreDiscrete {

	public KStarScoreUpperBound(MSKStarSettings settings) {
		super(settings);
	}

	public KStarScoreUpperBound(MSKStarSettings settings, PartitionFunction[] pfs) {
		super(settings, pfs);
	}

	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function lower bounds, so check 
		//against state-specific constraints that are upper bounds
		if(state <= numStates-2) {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, false);
		}

		//bound state partition function is an upper bound, so check 
		//against state-specific constraints that are lower bounds
		else {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, true);
		}
	}

	@Override
	public void compute(int maxNumConfs) {
		super.compute(maxNumConfs);
		/*
		//complex is null only for root node
		PartitionFunction complex = partitionFunctions[numStates-1];

		//here, the numerator is an upper bound, which only decreases
		//as we descend the tree. therefore, if numerator is 0, then we
		//can safely prune this node.
		if(complex != null && complex.getValues().qstar.compareTo(BigDecimal.ZERO)==0)
			constrSatisfied = false;
		 */
	}

	public BigDecimal getScore() {
		PartitionFunction complex = partitionFunctions[numStates-1];
		if(complex == null) return KStarScore.MAX_VALUE;

		BigDecimal score = super.getScore();

		//denom > 0
		if(getDenom().compareTo(BigDecimal.ZERO)>0) {
			return score;
		}

		//denom = 0
		else {
			//denom = 0. here, denom consists of lower bounds, which increase
			//as we descend the tree. so descendants might have good k* scores,
			//so we want to expand this node.
			return KStarScore.MAX_VALUE;
		}
	}

	public BigDecimal getUpperBoundScore() {
		return getScore();
	}

	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

}
