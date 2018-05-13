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

import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class BinaryRCSplitter extends RCSplitter {
	
	@Override
	public List<RC> split(int pos, RC rc) {
		
		// get the widest dof
		int dofi = -1;
		double maxWidth = 0;
		for (int i=0; i<rc.DOFs.size(); i++) {
			double width = rc.DOFmax.get(i) - rc.DOFmin.get(i);
			if (width > maxWidth) {
				maxWidth = width;
				dofi = i;
			}
		}
		
		double min = rc.DOFmin.get(dofi);
		double max = rc.DOFmax.get(dofi);
		double mid = (min + max)/2;
		
		return Arrays.asList(
			makeRC(rc, dofi, min, mid),
			makeRC(rc, dofi, mid, max)
		);
	}
}
