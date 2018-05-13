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

package edu.duke.cs.osprey.astar.conf.order;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.tools.MathTools;

public class DynamicHMeanAStarOrder implements AStarOrder {

	public final MathTools.Optimizer optimizer;

	public DynamicHMeanAStarOrder() {
		this(MathTools.Optimizer.Minimize);
	}

	public DynamicHMeanAStarOrder(MathTools.Optimizer optimizer) {
		this.optimizer = optimizer;
	}

	private AStarScorer gscorer;
	private AStarScorer hscorer;
	
	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		
		int bestPos = -1;
		double bestScore = optimizer.initDouble();

		for (int i=0; i<confIndex.numUndefined; i++) {
			
			int pos = confIndex.undefinedPos[i];
			double score = scorePos(confIndex, rcs, pos);

			if (optimizer.isBetter(score, bestScore)) {
				bestScore = score;
				bestPos = pos;
			}
		}

		if (bestPos >= 0) {
			return bestPos;
		}

		// sometimes, all the positions have infinite energies
		// so just pick one arbitrarily
		return confIndex.undefinedPos[0];
	}

	double scorePos(ConfIndex confIndex, RCs rcs, int pos) {
		
		// check all the RCs at this pos and aggregate the energies
		double parentScore = confIndex.node.getScore();
		double reciprocalSum = 0;
		for (int rc : rcs.get(pos)) {
			double childScore = gscorer.calcDifferential(confIndex, rcs, pos, rc)
				+ hscorer.calcDifferential(confIndex, rcs, pos, rc);
			reciprocalSum += 1.0/(childScore - parentScore);
		}

		// negate scores so better scores are lower, like energies
		return -1.0/reciprocalSum;
	}
}
