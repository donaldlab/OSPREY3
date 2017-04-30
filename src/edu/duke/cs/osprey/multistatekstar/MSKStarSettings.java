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
		
		//must deep copy, since these values can change
		this.search = new MSSearchProblem[other.search.length];
		System.arraycopy(other.search, 0, this.search, 0, other.search.length);
		
		this.isFinal = other.isFinal;
		this.constraints = other.constraints;
		
		this.pfTypes = new PartitionFunctionType[other.pfTypes.length];
		System.arraycopy(other.pfTypes, 0, this.pfTypes, 0, other.pfTypes.length);		
		
		this.ecalcs = new ConfEnergyCalculator.Async[other.ecalcs.length];
		System.arraycopy(other.ecalcs, 0, this.ecalcs, 0, other.ecalcs.length);
	}
}
