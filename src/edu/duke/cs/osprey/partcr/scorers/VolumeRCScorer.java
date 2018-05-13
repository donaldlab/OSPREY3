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

package edu.duke.cs.osprey.partcr.scorers;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.partcr.SplitWorld;

public class VolumeRCScorer implements RCScorer {

	@Override
	public double calcScore(SplitWorld splitWorld, RC rc, double boundErr) {
		
		// compare the voxel volume to the parent voxel volume
		RC parentRc = splitWorld.getSplits().getRCInfo(rc).getParent();
		assert (rc.DOFs.size() == parentRc.DOFs.size());
		int numDofs = rc.DOFs.size();
		
		double childVolume = 1;
		double parentVolume = 1;
		for (int i=0; i<numDofs; i++) {
			childVolume *= rc.DOFmax.get(i) - rc.DOFmin.get(i);
			parentVolume *= parentRc.DOFmax.get(i) - parentRc.DOFmin.get(i);
		}
		
		return boundErr*childVolume/parentVolume;
	}
}
