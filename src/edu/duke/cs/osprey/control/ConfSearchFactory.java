package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public interface ConfSearchFactory {
	
	ConfSearch make(EnergyMatrix emat, PruningMatrix pmat);
}
