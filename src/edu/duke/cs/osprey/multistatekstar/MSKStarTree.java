package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;

import edu.duke.cs.osprey.astar.AStarTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.confspace.SearchProblem;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class MSKStarTree extends AStarTree<FullAStarNode> {

	int numTreeLevels;//number of residues with sequence
	//changes+1 level if we are doing continuous minimization

	LMV objFcn;//we are minimizing objFcn
	LMV[] constraints;
	LMV[][] stateConstraints;

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;
	// MultiStateKStarTreeNode.assignments Assigns each level an index in 
	// AATypeOptions.get(level), and thus an AA type
	//If -1, then no assignment yet

	int numMaxMut;//number of mutations allowed away from wtSeq (-1 means no cap)
	ArrayList<String[]> wtSeqs;//bound state wild type sequences for each state

	int numStates;//how many states there are
	//states have the same mutable residues & options for AA residues,
	//but not necessarily for non AA residues

	SearchProblem search[][];//SearchProblems describing them
	//each state has >= 3 SearchProblems

	ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums;
	//mutable2StatePosNum.get(state) maps levels in this tree to flexible 
	//positions for state (not necessarily an onto mapping)

	int numSeqsReturned;

	public MSKStarTree(
			int numTreeLevels, 
			LMV objFcn, 
			LMV[] constraints,
			ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions, 
			int numMaxMut, 
			ArrayList<String[]> wtSeqs, 
			int numStates,
			SearchProblem[][] search, 
			ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums, 
			int numTopConfs) {

		this.numTreeLevels = numTreeLevels;
		this.objFcn = objFcn;
		this.constraints = constraints;
		this.AATypeOptions = AATypeOptions;
		this.numMaxMut = numMaxMut;
		this.wtSeqs = wtSeqs;
		this.numStates = numStates;
		this.search = search;
		this.mutable2StateResNums = mutable2StateResNums;
		
		numSeqsReturned = 0;
	}

	@Override
	public BigInteger getNumConformations() {
		throw new UnsupportedOperationException();
	}

	@Override
	public ArrayList<FullAStarNode> getChildren(FullAStarNode curNode) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public FullAStarNode rootNode() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isFullyAssigned(FullAStarNode node) {
		// TODO Auto-generated method stub
		return false;
	}
	
	public String seqAsString(int[] seqNodeAssignments) {
		return null;
	}

}
