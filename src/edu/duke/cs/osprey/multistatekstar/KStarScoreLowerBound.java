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

public class KStarScoreLowerBound extends KStarScoreDiscrete {

	public KStarScoreLowerBound(MSKStarSettings settings) {
		super(settings);
	}
	
	@Override
	protected void compute(int state, int maxNumConfs) {
		super.compute(state, maxNumConfs);
		//all unbound states are partition function upper bounds, so check 
		//against state-specific constraints that are lower bounds
		if(state <= numStates-2) {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, true);
		}
		
		//bound state is a partition function lower bound, so check
		//against state-specific constraints that are upper bounds
		else {
			if(constrSatisfied)
				constrSatisfied = checkConstraints(state, false);
		}
	}

}
