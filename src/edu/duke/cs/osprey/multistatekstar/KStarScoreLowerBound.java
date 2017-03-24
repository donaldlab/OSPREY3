package edu.duke.cs.osprey.multistatekstar;

public class KStarScoreLowerBound extends KStarScoreDiscrete {

	public KStarScoreLowerBound(MSKStarSettings settings) {
		super(settings);
	}
	
	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function upper bounds, so check 
		//against state-specific constraints that are lower bounds
		if(state <= partitionFunctions.length-2) {
			if(constrSatisfied) {
				constrSatisfied = checkConstraints(state, true);
			}
		}
	}

}
