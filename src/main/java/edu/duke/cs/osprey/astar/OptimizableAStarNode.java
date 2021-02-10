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

package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

public interface OptimizableAStarNode {

	double getGScore();
	void setGScore(double val);
	double getHScore();
	void setHScore(double val);
	int getLevel();

	default double getScore() {
		return getGScore() + getHScore();
	}

	default double getGScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getGScore(), optimizer);
	}
	default void setGScore(double val, MathTools.Optimizer optimizer) {
		setGScore(Tools.optimizeScore(val, optimizer));
	}

	default double getHScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getHScore(), optimizer);
	}
	default void setHScore(double val, MathTools.Optimizer optimizer) {
		setHScore(Tools.optimizeScore(val, optimizer));
	}

	default double getScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getScore(), optimizer);
	}

	public static class Tools {

		public static double optimizeScore(double score, MathTools.Optimizer optimizer) {
			switch (optimizer) {
				case Minimize: return score; // the pq is naturally a min-heap
				case Maximize: return -score; // negate the score so the pq acts like a max-heap
				default: throw new UnpossibleError();
			}
		}
	}
}
