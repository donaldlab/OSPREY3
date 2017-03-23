package edu.duke.cs.osprey.multistatekstar;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	MSSearchProblem[][] searchCont;//upper bound search problem for states
	MSSearchProblem[][] searchDisc;//lower bound search problem for state
	final MSKStarTree tree;//has all required objects
	double score;
	
	public MSKStarNode(
			MSKStarTree tree,
			MSSearchProblem[][] searchCont, 
			MSSearchProblem[][] searchDisc
			) {
		this.searchCont = searchCont;
		this.searchDisc = searchDisc;
		this.tree = tree;
		score = 0;
	}
	
	public void getScore() {
		
	}
	
	public void setScore(double val) {
		score = val;
	}
}

