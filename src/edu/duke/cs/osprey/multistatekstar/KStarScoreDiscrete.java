package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class KStarScoreDiscrete extends KStarScoreMinimized {

	public KStarScoreDiscrete(MSKStarSettings settings) {
		super(settings);
	}

	@Override
	public BigDecimal getLowerBoundScore() {
		return getScore();
	}

	@Override
	public BigDecimal getUpperBoundScore() {
		return getScore();
	}
	
	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		
		//multiply q* by number of undefined confs
		if(!settings.search[state].isFullyAssigned()) {
			PartitionFunctionDiscrete pf = (PartitionFunctionDiscrete) partitionFunctions[state];
			pf.getValues().qstar = pf.getValues().qstar.multiply(numUndefinedConfs(state));
		}
	}
	
	/**
	 * Compute either the minimum or maximum number of conformations of any
	 * possible undefined sub-sequence 
	 */
	private BigDecimal numUndefinedConfs(int state) {
		BigDecimal ans = BigDecimal.ONE;
		
		MSSearchProblem search = settings.search[state];
		boolean minConfs = search.settings.energyLBs ? false : true;
		PruningMatrix pmat = search.pruneMat;
		
		for(int pos : search.getPosNums(false)) {
			
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
			ans = ans.multiply(BigDecimal.valueOf(unPrunedConfs+prunedConfs));
			
		}
		
		return ans;
	}

}
