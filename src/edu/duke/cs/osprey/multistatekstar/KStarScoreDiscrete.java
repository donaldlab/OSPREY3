/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

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
