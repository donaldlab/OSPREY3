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

package edu.duke.cs.osprey.astar.conf.linked;

import edu.duke.cs.osprey.astar.conf.ConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.externalMemory.Queue;

public class LinkedConfAStarFactory implements ConfAStarFactory {
	
	@Override
	public Queue<ConfAStarNode> makeQueue(RCs rcs) {
		return Queue.PriorityFactory.of(null);
	}
	
	@Override
	public LinkedConfAStarNode makeRootNode(int numPos) {
		return new LinkedConfAStarNode();
	}
}
