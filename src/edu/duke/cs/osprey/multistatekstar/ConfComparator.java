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

import java.util.Comparator;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

/**
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */

public class ConfComparator implements Comparator<ScoredConf> {

	@Override
	public int compare(ScoredConf o1, ScoredConf o2) {
		double e1 = o1 instanceof EnergiedConf ? ((EnergiedConf)o1).getEnergy() : o1.getScore();
		double e2 = o2 instanceof EnergiedConf ? ((EnergiedConf)o2).getEnergy() : o2.getScore();
		return e1 <= e2 ? 1 : -1;
	}

}
