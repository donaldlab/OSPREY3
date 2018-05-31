package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public class KStarScoreLowerBound extends KStarScoreDiscrete {

	public KStarScoreLowerBound(MSKStarSettings settings) {
		super(settings);
	}

	public KStarScoreLowerBound(MSKStarSettings settings, PartitionFunction[] pfs) {
		super(settings, pfs);
	}

	@Override
	protected void compute(int state, long maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function upper bounds, so check 
		//against state-specific constraints that are lower bounds
		if(state <= numStates-2) {			
			if(constrSatisfied) {
				constrSatisfied = checkConstraints(state, true);
			}
		}

		//bound state is a partition function lower bound, so check
		//against state-specific constraints that are upper bounds
		else {
			if(constrSatisfied) {
				constrSatisfied = checkConstraints(state, false);
			}
		}
	}

	@Override
	public void compute(long maxNumConfs) {
		super.compute(maxNumConfs);

		if(getDenominator().compareTo(BigDecimal.ZERO)==0) {
			//denom = 0. here, denom consists of upper bounds, which only get
			//smaller as we descend the tree, meaning that all subsequent 
			//upper bounds will be 0. So we can safely prune this node.
			constrSatisfied = false;
		}
	}

	public BigDecimal getScore() {
		PartitionFunction complex = partitionFunctions[numStates-1];
		if(complex == null) return toLog10(BigDecimal.ZERO);

		//denom = 0
		//this node must be pruned
		if(getDenominator().compareTo(BigDecimal.ZERO)==0) {
			if(constrSatisfied) throw new RuntimeException("ERROR: this node must be pruned");
			return toLog10(KStarScore.MAX_VALUE);
		}

		//denom > 0
		else {
			BigDecimal score = super.getScore();
			return score;
		}
	}

	public BigDecimal getUpperBoundScore() {
		return getScore();
	}

	public BigDecimal getLowerBoundScore() {
		return getScore();
	}
}
