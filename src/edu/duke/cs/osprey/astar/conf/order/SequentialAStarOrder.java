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

public class SequentialAStarOrder implements AStarOrder {
	
	// NOTE: there's basically no reason to ever use this except to benchmark better methods
	// use a more intelligent static ordering instead

	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		// don't care...
	}
	
	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {
		
		// easy peasy
		// eg, the root node has level 0, so expand pos 0 next
		return confIndex.node.getLevel();
	}
}
