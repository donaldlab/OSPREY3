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
