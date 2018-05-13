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

package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;

public class KStarScoreUpperBound extends KStarScoreLowerBound {

	public KStarScoreUpperBound(MSKStarSettings settings) {
		super(settings);
	}

	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function lower bounds, so check 
		//against state-specific constraints that are upper bounds
		if(state <= numStates-2) {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, false);
		}

		//bound state partition function is an upper bound, so check 
		//against state-specific constraints that are lower bounds
		else {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, true);
		}
	}

	public BigDecimal getScore() {
		BigDecimal score = super.getScore();
		if(score.compareTo(BigDecimal.ZERO)>0) return score;
		else {
			if(getDenom().compareTo(BigDecimal.ZERO)>0) return score;
			else {
				if(initialized[numStates-1]) return score; //upper bound partition function is also 0
				else return PartitionFunctionMinimized.MAX_VALUE;
			}
		}
	}
	
	public BigDecimal getUpperBoundScore() {
		return getScore();
	}

}
