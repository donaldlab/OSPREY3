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

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class SmartRCSplitter extends RCSplitter {
	
	// NOTE: in quick (and certinaly not exhaustive tests),
	// this splitter didn't perform as well as the Binary or NAry splitters
	
	private double[] dofBoundValues;
	private int baseDofIndex;

	public SmartRCSplitter(double[] dofBoundValues, int baseDofIndex) {
		this.dofBoundValues = dofBoundValues;
		this.baseDofIndex = baseDofIndex;
	}

	@Override
	public List<RC> split(int pos, RC rc) {
		
		List<RC> rcs = null;
		
		// split along the widest dof, separating the two points down the middle
		int widestDofIndex = getWidestDofIndex(rc);
		if (widestDofIndex >= 0) {
			
			double min = rc.DOFmin.get(widestDofIndex);
			double max = rc.DOFmax.get(widestDofIndex);
			double boundVal = getBoundVal(widestDofIndex);
			double currentVal = getCurrentVal(rc, widestDofIndex);
			double mid = (boundVal + currentVal)/2;
			
			assert (boundVal >= min && boundVal <= max);
			assert (currentVal >= min && currentVal <= max);
			
			rcs = new ArrayList<>();
			rcs.add(makeRC(rc, widestDofIndex, min, mid));
			rcs.add(makeRC(rc, widestDofIndex, mid, max));
		}
		
		return rcs;
	}
	
	private int getWidestDofIndex(RC rc) {
		
		double bestDist = 0;
		int bestIndex = -1;
		
		for (int i=0; i<rc.DOFs.size(); i++) {
			
			double dist = Math.abs(getBoundVal(i) - getCurrentVal(rc, i));
			
			if (dist > bestDist) {
				bestDist = dist;
				bestIndex = i;
			}
		}
		
		return bestIndex;
	}
	
	private double getBoundVal(int dofIndex) {
		return dofBoundValues[baseDofIndex + dofIndex];
	}
	
	private double getCurrentVal(RC rc, int dofIndex) {
		return rc.DOFs.get(dofIndex).getCurVal();
	}
}
