package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
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
	
	public MSKStarSettings(MSKStarSettings other) {
		this.isReportingProgress = other.isReportingProgress;
		this.targetEpsilon = other.targetEpsilon;
		this.state = other.state;
		this.numTopConfsToSave = other.numTopConfsToSave;
		this.cfp = other.cfp;
		this.scoreType = other.scoreType;
		this.search = other.search;
		this.isFinal = other.isFinal;
		this.constraints = other.constraints;
		this.pfTypes = other.pfTypes;
		this.ecalcs = other.ecalcs;
	}
}
