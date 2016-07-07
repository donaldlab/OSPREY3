package edu.duke.cs.osprey.partcr.scorers;

import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.partcr.SplitWorld;

public class SplitsRCScorer implements RCScorer {
	
	private double splitPenalty;
	
	public SplitsRCScorer() {
		this(1);
	}
	
	public SplitsRCScorer(double splitPenalty) {
		this.splitPenalty = splitPenalty;
	}
	
	@Override
	public double calcScore(SplitWorld splitWorld, RC rc, double boundErr) {
		
		int numVoxels = splitWorld.getSplits().getRCInfo(rc).getNumVoxels();
		return boundErr/(splitPenalty*(numVoxels - 1) + 1);
	}
}
