package edu.duke.cs.osprey.partcr.splitters;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class NAryRCSplitter extends RCSplitter {
	
	private static enum Side {
		
		Low,
		High;
		
		public static Side get(int rcIndex, int dofIndex) {
			// whats the nth bit set to?
			int ord = (rcIndex >> dofIndex) & 0x1;
			return values()[ord];
		}
	}
	
	@Override
	public List<RC> split(int pos, RC rc) {
		
		List<RC> rcs = new ArrayList<>();
		
		// get the splittable DoFs
		List<Integer> indices = new ArrayList<>();
		for (int i=0; i<rc.DOFs.size(); i++) {
			if (rc.DOFmin.get(i) != rc.DOFmax.get(i)) {
				indices.add(i);
			}
		}
		int numDofs = indices.size();
		
		// split all the dimensions at once
		int numRcs = 1 << numDofs;
		for (int rcIndex=0; rcIndex<numRcs; rcIndex++) {
			
			// clone the RC
			ArrayList<Double> mins = new ArrayList<>(rc.DOFmin);
			ArrayList<Double> maxs = new ArrayList<>(rc.DOFmax);
			RC subRc = new RC(rc.AAType, rc.template, rc.rotNum, rc.DOFs, mins, maxs, rc.RCIndex);
			
			// change the dof bounds
			for (int dofIndex : indices) {
				
				// split the dof down the middle
				double min = rc.DOFmin.get(dofIndex);
				double max = rc.DOFmax.get(dofIndex);
				double mid = (max + min)/2;
				
				switch (Side.get(rcIndex, dofIndex)) {
					case Low:
						mins.set(dofIndex, min);
						maxs.set(dofIndex, mid);
					break;
					case High:
						mins.set(dofIndex, mid);
						maxs.set(dofIndex, max);
					break;
				}
			}
			
			rcs.add(subRc);
		}
		
		// TEMP: check the split RCs
		double volume = calcVolume(rc, indices);
		for (RC subRc : rcs) {
			volume -= calcVolume(subRc, indices);
		}
		assert (Math.abs(volume) < 1e-10);
		
		return rcs;
	}
	
	private double calcVolume(RC rc, List<Integer> indices) {
		double volume = 1;
		for (int dofIndex : indices) {
			double min = rc.DOFmin.get(dofIndex);
			double max = rc.DOFmax.get(dofIndex);
			volume *= max - min;
		}
		return volume;
	}
}
