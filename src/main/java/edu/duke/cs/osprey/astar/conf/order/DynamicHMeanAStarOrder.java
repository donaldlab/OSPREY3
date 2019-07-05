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
	public boolean isDynamic() {
		return true;
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
