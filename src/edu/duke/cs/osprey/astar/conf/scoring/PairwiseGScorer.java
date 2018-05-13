/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.tools.MathTools;

public class PairwiseGScorer implements AStarScorer {
	
	public final EnergyMatrix emat;
	public final MathTools.Optimizer optimizer;

	public PairwiseGScorer(EnergyMatrix emat) {
		this(emat, MathTools.Optimizer.Minimize);
	}

	public PairwiseGScorer(EnergyMatrix emat, MathTools.Optimizer optimizer) {
		this.emat = emat;
		this.optimizer = optimizer;
	}
	
	@Override
	public PairwiseGScorer make() {
		return new PairwiseGScorer(emat, optimizer);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// constant term
    	double gscore = emat.getConstTerm();
    	
    	// one body energies
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];
			
			gscore += emat.getOneBody(pos1, rc1);
		}
		
		// pairwise energies
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];
		
			for (int j=0; j<i; j++) {
				int pos2 = confIndex.definedPos[j];
				int rc2 = confIndex.definedRCs[j];
				
				gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
		return gscore;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
    	// modify the parent node's g-score
    	double gscore = confIndex.node.getGScore(optimizer);
    	
    	// add the new one-body energy
    	gscore += emat.getOneBody(nextPos, nextRc);
    	
    	// add the new pairwise energies
    	for (int i=0; i<confIndex.numDefined; i++) {
    		int pos = confIndex.definedPos[i];
    		int rc = confIndex.definedRCs[i];
    		gscore += emat.getPairwise(pos, rc, nextPos, nextRc);
    	}
    	
    	return gscore;
	}
}
