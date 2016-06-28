package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.tools.HashCalculator;

public class NAryRCSplitter implements RCSplitter {
	
	private static enum Side {
		
		Low,
		High;
		
		public static Side get(int rcIndex, int dofIndex) {
			// whats the nth bit set to?
			int ord = (rcIndex >> dofIndex) & 0x1;
			return values()[ord];
		}
	}
	
	private static class RotId {
		
		private int pos;
		private int rotNum;
		
		public RotId(int pos, int rotNum) {
			this.pos = pos;
			this.rotNum = rotNum;
		}
		
		@Override
		public boolean equals(Object other) {
			if (other instanceof RotId) {
				return equals((RotId)other);
			}
			return false;
		}
		
		public boolean equals(RotId other) {
			return this.pos == other.pos
				&& this.rotNum == other.rotNum;
		}
		
		@Override
		public int hashCode() {
			return HashCalculator.hashIds(pos, rotNum);
		}
	}
	
	
	private Set<RotId> splitRotIds;
	
	public NAryRCSplitter() {
		splitRotIds = new HashSet<>();
	}
	
	@Override
	public boolean willSplit(int pos, RC rc) {
		
		// only split each rotamer once
		return !splitRotIds.contains(new RotId(pos, rc.rotNum));
	}

	@Override
	public List<RC> split(int pos, RC rc) {
		
		// make sure we don't split this one again
		splitRotIds.add(new RotId(pos, rc.rotNum));
		
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
