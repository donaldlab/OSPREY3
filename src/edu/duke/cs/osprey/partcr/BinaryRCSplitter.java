package edu.duke.cs.osprey.partcr;

import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.confspace.RC;

public class BinaryRCSplitter extends RCSplitter {
	
	private double minWidth;
	
	public BinaryRCSplitter() {
		this(0);
	}
	
	public BinaryRCSplitter(double minWidth) {
		this.minWidth = minWidth;
	}

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
		
		// don't keep splitting if we hit the floor
		if (maxWidth <= minWidth) {
			return null;
		}
		
		double min = rc.DOFmin.get(dofi);
		double max = rc.DOFmax.get(dofi);
		double mid = (min + max)/2;
		
		// TEMP
		System.out.println(String.format("split RC %d-%d along dof %d: [%f,%f] %f",
			pos, rc.rotNum, dofi, min, max, mid
		));
		
		return Arrays.asList(
			makeRC(rc, dofi, min, mid),
			makeRC(rc, dofi, mid, max)
		);
	}
}
