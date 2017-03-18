package edu.duke.cs.osprey.multistatekstar;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	MSSearchProblem[] discSearch;//lower bound search problem for state
	MSSearchProblem[] contSearch;//upper bound search problem for state
	MSKStarScore[] lbScores;//lower bound k* scores for state
	MSKStarScore[] ubScores;//upper bound k* scores for state
	
	public MSKStarNode(MSSearchProblem[] discSearch, MSSearchProblem[] contSearch) {
		this.discSearch = discSearch;
		this.contSearch = contSearch;
	}
	
	/**
	 * compute upper and lower bound K* scores.
	 * compute node score.
	 */
	public void getScore() {
		
	}
}

