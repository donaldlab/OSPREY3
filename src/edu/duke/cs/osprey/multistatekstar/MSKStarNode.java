package edu.duke.cs.osprey.multistatekstar;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	MSSearchProblem[] lbSearch;//lower bound pruning matrix for each state
	MSSearchProblem[] ubSearch;//upper bound pruning matrix for each state
	KStarScore[] lbScores;//lower bound k* scores for each state
	KStarScore[] ubScores;//upper bound k* scores for each state
	
	public MSKStarNode(MSSearchProblem[] lbSearch, MSSearchProblem[] ubSearch) {
		this.lbSearch = lbSearch;
		this.ubSearch = ubSearch;
	}
}
