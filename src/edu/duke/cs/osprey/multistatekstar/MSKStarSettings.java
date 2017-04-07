package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.KStarScore.PartitionFunctionType;
/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public class MSKStarSettings {
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public int state;
	public int numTopConfsToSave;
	public MSConfigFileParser cfp;
	public KStarScoreType scoreType;
	public MSSearchProblem[] search;
	public boolean isFinal;
	public LMB[] constraints;
	public PartitionFunctionType[] pfTypes;
	public ConfEnergyCalculator.Async[] ecalcs;

	public MSKStarSettings() {}
}
