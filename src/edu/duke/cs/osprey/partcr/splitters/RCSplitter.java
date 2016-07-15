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
