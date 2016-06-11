package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public interface MPLPUpdater {
	
	void update(MessageVars lambdas, EnergyMatrix emat); 
}
