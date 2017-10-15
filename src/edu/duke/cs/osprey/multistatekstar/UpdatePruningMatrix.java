package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.pruning.PruningMatrix;

public interface UpdatePruningMatrix {

	public PruningMatrix updatePruningMatrix( ArrayList<Integer> splitPosNums, ArrayList<ArrayList<String>> splitAAs);
	
}
