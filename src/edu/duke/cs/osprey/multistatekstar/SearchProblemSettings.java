package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class SearchProblemSettings {

	public PruningMatrix invPruneMat;//inverted pruning matrix
	public PruningMatrix redPruneMat;//pruning matrix reduced to new posnums and allowed AAs
	
	public EnergyMatrix redEmat;//emat reduced to new posnums and allowed AAs

	public ArrayList<String> redFlexRes;//reduced flex res
	public ArrayList<ArrayList<String>> redAllowedAAs;//reduced allowed AAs
	
	public ParamSet params;
	
	public SearchProblemSettings() {}
}
