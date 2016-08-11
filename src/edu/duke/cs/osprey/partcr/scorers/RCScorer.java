package edu.duke.cs.osprey.partcr.scorers;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.partcr.SplitWorld;

public interface RCScorer {
	
	double calcScore(SplitWorld splitWorld, RC rc, double boundErr);
}
