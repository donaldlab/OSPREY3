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
