package edu.duke.cs.osprey.multistatekstar;

public class KStarScoreUpperBound extends KStarScoreLowerBound {

	public KStarScoreUpperBound(MSKStarSettings settings) {
		super(settings);
	}

	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function lower bounds, so check 
		//against state-specific constraints that are upper bounds
		if(state <= partitionFunctions.length-2) {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, false);
		}

		//bound state partition function is an upper bound, so check 
		//against state-specific constraints that are lower bounds
		else if(state == partitionFunctions.length-1) {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, true);
		}
	}

}
