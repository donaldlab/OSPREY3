package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class SmartRCSplitter extends RCSplitter {
	
	private double[] dofBoundValues;
	private int baseDofIndex;

	public SmartRCSplitter(double[] dofBoundValues) {
		this.dofBoundValues = dofBoundValues;
		this.baseDofIndex = 0;
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
			
			// TEMP
			System.out.println(String.format("split RC along dof %d: min:%f, bound:%f, [%f,%f] %f", widestDofIndex, currentVal, boundVal, min, max, mid));
			assert (boundVal >= min && boundVal <= max);
			assert (currentVal >= min && currentVal <= max);
			
			rcs = new ArrayList<>();
			rcs.add(makeRC(rc, widestDofIndex, min, mid));
			rcs.add(makeRC(rc, widestDofIndex, mid, max));
		}
		
		baseDofIndex += rc.DOFs.size();
		
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
