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

package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.TuplesIndex;


public class LUTEGScorer implements AStarScorer {

	public final LUTEConfEnergyCalculator ecalc;

	public LUTEGScorer(LUTEConfEnergyCalculator ecalc) {
		this.ecalc = ecalc;
	}

	@Override
	public AStarScorer make() {
		return new LUTEGScorer(ecalc);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {

		/* convert the conf index into a conf
			yeah, we'll take a slight performance hit doing this,
			but the higher-order tuple indices are optimized around fast pos->rc lookups,
			but ConfIndex optimizes fast singles/pairs enumeration instead
		*/
		int[] conf = Conf.make(confIndex);

		try {
			return ecalc.calcEnergy(conf);
		} catch (TuplesIndex.NoSuchTupleException ex) {
			// conf has a pruned tuple, can't score it
			return Double.POSITIVE_INFINITY;
		}
	}
}
