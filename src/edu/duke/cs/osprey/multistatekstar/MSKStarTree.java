package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.multistatekstar.ResidueOrder.ResidueOrderType;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class MSKStarTree {

	int numTreeLevels;//number of residues with sequence
	//changes+1 level if we are doing continuous minimization

	LMV objFcn;//we are minimizing objFcn
	LMV[] constraints;
	LMV[][] stateConstraints;

	ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums;
	//mutable2StateResNum.get(state) maps levels in this tree to flexible positions for state

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;
	// MultiStateKStarTreeNode.assignments Assigns each level an index in 
	// AATypeOptions.get(level), and thus an AA type
	//If -1, then no assignment yet

	int numMaxMut;//number of mutations allowed away from wtSeq (-1 means no cap)
	ArrayList<String[]> wtSeqs;//bound state wild type sequences for each state

	int numStates;//how many states there are
	//states have the same mutable residues & options for AA residues,
	//but not necessarily for non AA residues

	SearchProblem searchCont[][];//SearchProblems describing them; each state has >= 3 SearchProblems
	SearchProblem searchDisc[][];

	ConfEnergyCalculator.Async[][] ecalcsCont;//energy calculators for continuous emats
	ConfEnergyCalculator.Async[][] ecalcsDisc;//energy calculators for discrete emats

	ParamSet msParams;//multistate spec params
	MSConfigFileParser[] cfps;//config file parsers for each state

	ResidueOrderType residueOrder;

	int numSeqsReturned;

	public MSKStarTree(
			int numTreeLevels,
			int numStates,
			int numMaxMut,
			LMV objFcn,
			LMV[] constraints,
			ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums,
			ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions,  
			ArrayList<String[]> wtSeqs, 
			SearchProblem[][] searchCont,
			SearchProblem[][] searchDisc,
			ConfEnergyCalculator.Async[][] ecalcsCont,
			ConfEnergyCalculator.Async[][] ecalcsDisc,
			ParamSet msParams,
			MSConfigFileParser[] cfps
			) {

		this.numTreeLevels = numTreeLevels;
		this.objFcn = objFcn;
		this.constraints = constraints;
		this.AATypeOptions = AATypeOptions;
		this.numMaxMut = numMaxMut;
		this.wtSeqs = wtSeqs;
		this.numStates = numStates;
		this.searchCont = searchCont;
		this.searchDisc = searchDisc;
		this.ecalcsCont = ecalcsCont;
		this.ecalcsDisc = ecalcsDisc;
		this.mutable2StateResNums = mutable2StateResNums;

		this.cfps = cfps;
		this.msParams = msParams;
		residueOrder = getResidueOrder();
		numSeqsReturned = 0;
	}

	private ResidueOrderType getResidueOrder() {
		String val = msParams.getValue("RESIDUEORDER");
		switch(val.toLowerCase()) {
		case "staticsequential":
			return ResidueOrderType.StaticSequential;
		case "staticmindom":
			return ResidueOrderType.StaticMinDom;
		case "staticobjFunchmean":
			return ResidueOrderType.StaticObjFuncHMean;
		case "dynamicbbjfunchmean":
			return ResidueOrderType.DynamicObjFuncHMean;
		default:
			throw new UnsupportedOperationException("ERROR: unsupported residue order type: "+val);
		}
	}

	public ArrayList<MSKStarNode> getChildren(MSKStarNode curNode) {
		//for each state and substate, pick next position to expand

		//create search problems

		//compute and score children
		ArrayList<MSKStarNode> ans = new ArrayList<>();
		ans.trimToSize();
		return ans;
	}

	public MSKStarNode rootNode() {
		// TODO Auto-generated method stub
		return null;
	}

	public boolean isLeafNode(MSKStarNode node) {
		// TODO Auto-generated method stub
		return false;
	}

	public String[] nextSeq() {
		//a sequence spans multiple states
		return null;
	}

}
