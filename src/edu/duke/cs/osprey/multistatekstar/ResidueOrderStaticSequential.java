package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class ResidueOrderStaticSequential extends ResidueOrder {

	public ResidueOrderStaticSequential() {
		super();
	}

	private ArrayList<ArrayList<AAAssignment>> getBoundAAAssignments(
			MSSearchProblem[] search, 
			int splitPos, 
			boolean getPosFromUnbound, 
			int numMaxMut) {

		ArrayList<Integer> complexPos = new ArrayList<>();
		MSSearchProblem complex = search[search.length-1];

		//root node, so the unbound states give us bound state splitpos
		if(getPosFromUnbound) {
			for(int subState=0;subState<search.length-1;++subState)
				complexPos.add(complex.flexRes.indexOf(search[subState].flexRes.get(splitPos)));
		}

		//not root node, so adding a single pos in the bound state
		else complexPos.add(splitPos);

		complexPos.trimToSize();

		ArrayList<ArrayList<AAAssignment>> ans = new ArrayList<>();
		String[] wt = MSKStarNode.WT_SEQS.get(0);
		String[] buf = new String[wt.length];
		getBoundAAAssignmentsHelper(complex.settings.AATypeOptions, ans, complexPos, wt, buf, 0, 0, numMaxMut);

		ans.trimToSize();
		return ans;
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextResidueAssignment(LMB objFcn,
			MSSearchProblem[][] objFcnSearch, int numMaxMut) {

		int state = 0;
		int numSubStates = objFcnSearch[state].length;
		//bound state is the sequence
		MSSearchProblem boundState = objFcnSearch[state][numSubStates-1];

		if(boundState.isFullyAssigned())//no positions to split
			throw new RuntimeException("ERROR: there are no unassigned positions");

		ArrayList<ArrayList<AAAssignment>> boundAssignments = null;
		if(boundState.getNumAssignedPos()==0) {//root node; add all allowed single mutations from unbound states
			boundAssignments = getBoundAAAssignments(objFcnSearch[state], 0, true, numMaxMut);
		}

		else {//add all allowed mutations at the next numerical splitPos
			ArrayList<Integer> splitPos = boundState.getPosNums(false);//sequential = first unassigned pos
			boundAssignments = getBoundAAAssignments(objFcnSearch[state], splitPos.get(0), false, numMaxMut);
		}

		ArrayList<ArrayList<ArrayList<AAAssignment>>> ans = new ArrayList<>();
		//splitting a bound state position means splitting the corresponding unbound state position(s)
		for(int subState=0;subState<numSubStates-1;++subState) {
			MSSearchProblem unbound = objFcnSearch[state][subState];
			ans.add(getUnboundAAAssignments(boundState, boundAssignments, unbound));
		}
		ans.add(boundAssignments);

		ans.trimToSize();
		return ans;
	}
}
