package edu.duke.cs.osprey.partcr.scorers;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.partcr.SplitWorld;

public class NopRCScorer implements RCScorer {

	@Override
	public double calcScore(SplitWorld splitWorld, RC rc, double boundErr) {
		
		// no-op, just return the error
		return boundErr;
	}
}
