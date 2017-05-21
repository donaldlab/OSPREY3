package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class KStarScoreDiscrete extends KStarScoreMinimized {

	public KStarScoreDiscrete(MSKStarSettings settings) {
		super(settings);
	}

	public KStarScoreDiscrete(MSKStarSettings settings, PartitionFunction[] pfs) {
		super(settings, pfs);
	}

	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		return super.getUpperBoundScore();
	}

	@Override
	public void compute(int maxNumConfs) {
		super.compute(maxNumConfs);
	}

	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
	
		PartitionFunctionDiscrete pf = (PartitionFunctionDiscrete) partitionFunctions[state];
		
		//multiply by assigned*unassigned confs
		if(settings.pfTypes[state] == PartitionFunctionType.UpperBound) {
			BigDecimal unassignedConfs = numConfs(state, false);
			BigDecimal assignedConfs = numConfs(state, true);
			pf.getValues().qstar = pf.getValues().qstar.multiply(assignedConfs.multiply(unassignedConfs));
		}
		
		//only multiply by unassigned confs
		else if(!settings.search[state].isFullyAssigned()) {	
			BigDecimal unassignedConfs = numConfs(state, false);
			pf.getValues().qstar = pf.getValues().qstar.multiply(unassignedConfs);
		}
	}

	/**
	 * Compute either the minimum or maximum number of conformations of any
	 * possible undefined sub-sequence 
	 */
	private BigDecimal numConfs(int state, boolean assigned) {
		BigDecimal ans = BigDecimal.ONE;

		MSSearchProblem search = settings.search[state];
		boolean minConfs = search.settings.energyLBs ? false : true;
		PruningMatrix pmat = search.pruneMat;

		for(int pos : search.getPosNums(assigned)) {

			long unPrunedConfs = minConfs ? Long.MAX_VALUE : Long.MIN_VALUE;
			long prunedConfs = minConfs ? Long.MAX_VALUE : Long.MIN_VALUE;

			for(String AAType : search.settings.AATypeOptions.get(pos)) {
				long numAARCs = search.unprunedAtPos(pmat, pos, AAType).size();
				unPrunedConfs = minConfs ? Math.min(unPrunedConfs, numAARCs) : Math.max(unPrunedConfs, numAARCs);
				if(!minConfs) {
					numAARCs = search.unprunedAtPos(partitionFunctions[state].invmat, pos, AAType).size();
					prunedConfs = Math.max(prunedConfs, numAARCs);
				}
			}

			if(minConfs) prunedConfs = 0;
			//corner case where all confs are pruned due to mutation's intrinsic steric clash
			if(!minConfs && (unPrunedConfs+prunedConfs)==0) continue;
			ans = ans.multiply(BigDecimal.valueOf(unPrunedConfs+prunedConfs));
		}

		return ans;
	}

}
