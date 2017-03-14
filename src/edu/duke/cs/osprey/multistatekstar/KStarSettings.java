package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.gmec.ConfEnergyCalculator;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.KStarScore.PartitionFunctionType;

public class KStarSettings {
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public int state;
	public int numTopConfsToSave;
	public MSConfigFileParser cfp;
	public KStarScoreType scoreType;
	public MSSearchProblem[] search;
	public LMV[] constraints;
	public PartitionFunctionType[] pfTypes;
	public ConfEnergyCalculator.Async[] ecalcs;

	public KStarSettings() {}
}
