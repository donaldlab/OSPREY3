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

package edu.duke.cs.osprey.astar.conf.pruning;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;

public interface AStarPruner {

	/**
	 * check to see if a node should be pruned, after it comes off the heap
	 */
	boolean isPruned(ConfAStarNode node);

	/**
	 * check to see if a child node should be pruned before it's created and scored
	 */
	boolean isPruned(ConfAStarNode node, int nextPos, int nextRc);
}
