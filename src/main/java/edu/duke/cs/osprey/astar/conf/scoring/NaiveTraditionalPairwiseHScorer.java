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

public class NaiveTraditionalPairwiseHScorer implements AStarScorer {
	
	private EnergyMatrix emat;
	
	public NaiveTraditionalPairwiseHScorer(EnergyMatrix emat) {
		this.emat = emat;
	}
	
	@Override
	public NaiveTraditionalPairwiseHScorer make() {
		return new NaiveTraditionalPairwiseHScorer(emat);
	}
	
	@Override
	public double calc(ConfIndex index, RCs rcs) {
		
		double hscore = 0;
		
		// get the score for each undefined position
		for (int i=0; i<index.numUndefined; i++) {
			int pos1 = index.undefinedPos[i];
			
			// min over possible assignments to pos1
			double pos1Score = Double.POSITIVE_INFINITY;
			for (int rc1 : rcs.get(pos1)) {
				
				double rcContrib = emat.getOneBody(pos1, rc1);

				// interactions with defined residues
				for (int j=0; j<index.numDefined; j++) {
					int pos2 = index.definedPos[j];
					int rc2 = index.definedRCs[j];
					rcContrib += emat.getPairwise(pos1, rc1, pos2, rc2);
				}

				// interactions with undefined residues
				for (int j=0; j<index.numUndefined; j++) {
					int pos2 = index.undefinedPos[j];
					if (pos2 >= pos1) {
						break;
					}

					// min over possible assignments to pos2
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
						minEnergy = Math.min(minEnergy, pairwiseEnergy);
					}

					rcContrib += minEnergy;
				}

				pos1Score = Math.min(pos1Score, rcContrib);
			}
		
			hscore += pos1Score;
		}
		
		return hscore;
	}
}
