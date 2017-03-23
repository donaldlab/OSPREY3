package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class DiscreteKStarScore extends ContinuousKStarScore {

	public DiscreteKStarScore(MSKStarSettings settings) {
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
		if(!settings.search[state].isFullyDefined()) {
			DiscretePartitionFunction pf = (DiscretePartitionFunction) partitionFunctions[state];
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
		boolean min = search.contSCFlex ? false : true;
		QPruningMatrix pmat = (QPruningMatrix)search.pruneMat;
		
		for(int pos : search.getPos(false)) {
			
			long unPrunedConfs = min ? Long.MAX_VALUE : Long.MIN_VALUE;
			long prunedConfs = min ? Long.MAX_VALUE : Long.MIN_VALUE;
			
			for(String AAType : search.allowedAAs.get(pos)) {
				long numAARCs = search.unprunedAtPos(pmat, pos, AAType).size();
				unPrunedConfs = min ? Math.min(unPrunedConfs, numAARCs) : Math.max(unPrunedConfs, numAARCs);
				if(!min) {
					numAARCs = search.unprunedAtPos((QPruningMatrix)pmat.invert(), pos, AAType).size();
					prunedConfs = Math.max(prunedConfs, numAARCs);
				}
			}
			
			if(min) prunedConfs = 0;
			ans = ans.multiply(BigDecimal.valueOf(unPrunedConfs+prunedConfs));
			
		}
		
		return ans;
	}

}
