package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public class KStarScoreUpperBound extends KStarScoreLowerBound {

	public KStarScoreUpperBound(MSKStarSettings settings) {
		super(settings);
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

	public BigDecimal getScore() {
		BigDecimal score = super.getScore();
		if(score.compareTo(BigDecimal.ZERO)>0) return score;
		else {
			if(getDenom().compareTo(BigDecimal.ZERO)>0) return score;
			else {
				if(initialized[numStates-1]) return score; //upper bound partition function is also 0
				else return PartitionFunctionMinimized.MAX_VALUE;
			}
		}
	}
	
	public BigDecimal getUpperBoundScore() {
		return getScore();
	}

}
