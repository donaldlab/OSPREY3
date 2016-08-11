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
