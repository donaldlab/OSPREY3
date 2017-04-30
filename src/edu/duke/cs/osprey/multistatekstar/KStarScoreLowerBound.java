package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public class KStarScoreLowerBound extends KStarScoreDiscrete {

	public KStarScoreLowerBound(MSKStarSettings settings) {
		super(settings);
	}
	
	public KStarScoreLowerBound(MSKStarSettings settings, PartitionFunction[] pfs) {
		super(settings, pfs);
	}

	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function upper bounds, so check 
		//against state-specific constraints that are lower bounds
		if(state <= numStates-2) {			
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, true);
		}
		
		//bound state is a partition function lower bound, so check
		//against state-specific constraints that are upper bounds
		else {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, false);
		}
	}

}
