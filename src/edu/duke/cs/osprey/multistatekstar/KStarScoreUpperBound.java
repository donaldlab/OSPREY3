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