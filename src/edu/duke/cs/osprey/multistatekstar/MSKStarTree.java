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
public class MSKStarTree extends AStarTree<FullAStarNode> {

	int numTreeLevels;//number of residues with sequence changes+1
	//level if we are doing continuous minimization

	LMV objFcn;//we are minimizing objFcn
	LMV[] kssConstraints;
	LMV[][] pfConstraints;

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;
	// MultiStateKStarTreeNode.assignments Assigns each level an index in 
	// AATypeOptions.get(level), and thus an AA type
	//If -1, then no assignment yet

	int numMaxMut;//number of mutations allowed away from wtSeq (-1 means no cap)
	String wtSeq[];//wild type sequence

	//information on states

	int numStates;//how many states there are
	//states have the same mutable residues & options for AA residues,
	//but not necessarily for non AA residues

	SearchProblem search[][];//SearchProblems describing them
	//each state has >= 3 SearchProblems

	ArrayList<ArrayList<Integer>> mutable2StatePosNums;
	//mutable2StatePosNum.get(state) maps levels in this tree to flexible 
	//positions for state (not necessarily an onto mapping)

	int stateNumPos[];

	int numSeqsReturned = 0;
	int stateGMECsForPruning = 0;//how many state GMECs have been calculated for nodes that are pruned


	boolean outputGMECStructs;//Output GMEC structure for each (state, sequence)

	private static final long serialVersionUID = 1L;

	public MSKStarTree(
			int numTreeLevels, 
			LMV objFcn, 
			LMV[] constraints,
			ArrayList<ArrayList<String>> AATypeOptions, 
			int numMaxMut, 
			String[] wtSeq, 
			int numStates,
			SearchProblem[] stateSP, 
			ArrayList<ArrayList<ArrayList<Integer>>> mutable2StatePosNums, 
			int numTopConfs) {

	}

	@Override
	public BigInteger getNumConformations() {
		// TODO Auto-generated method stub
		return null;
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
