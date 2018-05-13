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
