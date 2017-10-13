package edu.duke.cs.osprey.multistatekstar;

import java.util.HashMap;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.KStarScore.PartitionFunctionType;
/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
public class MSKStarSettings {
	
	public static double TIMEOUT_HRS = Double.MAX_VALUE;
	public static HashMap<Integer, HashMap<Integer, Boolean>> MEMOIZE_STATE_PFS;
	
	public boolean isReportingProgress;
	public double targetEpsilon;
	public int state;
	public int numTopConfsToSave;
	public MSConfigFileParser cfp;
	public KStarScoreType scoreType;
	public MSSearchProblem[] search;
	public boolean isFinal;
	public boolean computeGMEC;
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
		this.computeGMEC = other.computeGMEC;
		
		this.search = new MSSearchProblem[other.search.length];
		System.arraycopy(other.search, 0, this.search, 0, other.search.length);
		
		this.isFinal = other.isFinal;
		this.constraints = other.constraints;
		
		this.pfTypes = new PartitionFunctionType[other.pfTypes.length];
		System.arraycopy(other.pfTypes, 0, this.pfTypes, 0, other.pfTypes.length);		
		
		this.ecalcs = new ConfEnergyCalculator.Async[other.ecalcs.length];
		System.arraycopy(other.ecalcs, 0, this.ecalcs, 0, other.ecalcs.length);
	}
	
	static boolean memoizePFs(int state, int substate) {
		if(MEMOIZE_STATE_PFS == null || MEMOIZE_STATE_PFS.get(state) == null) return false;
		Boolean ans = MEMOIZE_STATE_PFS.get(state).get(substate);
		return ans == null ? false : ans;
	}
}
