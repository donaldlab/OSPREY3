package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.control.ConfEnergyCalculator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	MSSearchProblem[][] searchCont;//upper bound search problem for states
	MSSearchProblem[][] searchDisc;//lower bound search problem for state
	ConfEnergyCalculator.Async[][] ecalcsCont;//continuous energy calculators
	ConfEnergyCalculator.Async[][] ecalcsDisc;//continuous energy calculators
	
	public MSKStarNode(
			MSSearchProblem[][] searchCont, 
			MSSearchProblem[][] searchDisc,
			ConfEnergyCalculator.Async[][] ecalcsCont,
			ConfEnergyCalculator.Async[][] ecalcsDisc,
			LMV objFcn,
			LMV[] constraints
			) {
		this.searchCont = searchCont;
		this.searchDisc = searchDisc;
		this.ecalcsCont = ecalcsCont;
		this.ecalcsDisc = ecalcsDisc;
	}
	
	/**
	 * compute upper and lower bound K* scores.
	 * compute node score.
	 */
	public void getScore() {
		
	}
}

