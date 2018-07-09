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
