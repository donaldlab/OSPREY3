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

package edu.duke.cs.osprey.partcr.splitters;

import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public abstract class RCSplitter {

	public abstract List<RC> split(int pos, RC rc);
	
	protected RC makeRC(RC rc, int dofIndex, double min, double max) {
		RC subRc = new RC(rc);
		subRc.DOFmin.set(dofIndex, min);
		subRc.DOFmax.set(dofIndex, max);
		return subRc;
	}
}
