package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Domain size determined by the number of rotamers in the bound state
 */
@SuppressWarnings("serial")
public class ResidueOrderStaticMinDomain extends ResidueOrderDynamicScore {

	private HashMap<Integer, Integer> residueValues;//residue values for bound state
	//in each unassigned pos
	private boolean computeProduct;

	public ResidueOrderStaticMinDomain(MSSearchProblem[][] objFcnSearch, boolean computeProduct) {
		super(objFcnSearch, true);
		this.residueValues = computeResidueValues(objFcnSearch);
		this.computeProduct = computeProduct;
	}

	private HashMap<Integer, Integer> computeResidueValues(MSSearchProblem[][] objFcnSearch) {
		int complex = objFcnSearch[0].length-1;
		ArrayList<Integer> unassignedPos = objFcnSearch[0][complex].getPosNums(false);

		HashMap<Integer, Integer> ans = new HashMap<>();

		for(int pos : unassignedPos) {
			int sumRCs = 0;
			for(int state=0;state<objFcnSearch.length;++state) {
				sumRCs += objFcnSearch[state][complex].pruneMat.unprunedRCsAtPos(pos).size();
			}
			ans.put(pos, sumRCs);
		}

		return ans;
	}
	
	public int residue2Domain(int residue) {
		return residueValues.get(residue);
	}

	//score is product of domain sizes of bound states
	protected BigDecimal getBoundStateDomainProduct(ResidueAssignment assignment) {
		long prodRCs = 1;
		int complex = assignment.length()-1;
		for(int pos : assignment.get(complex)) {
			prodRCs *= residueValues.get(pos);
		}
		return BigDecimal.valueOf(prodRCs);
	}

	//score is sum of domain sizes of bound states
	protected BigDecimal getBoundStateDomainSum(ResidueAssignment assignment) {
		long sumRCs = 1;
		int complex = assignment.length()-1;
		for(int pos : assignment.get(complex)) {
			sumRCs += residueValues.get(pos);
		}
		return BigDecimal.valueOf(sumRCs);
	}

	protected ArrayList<ResidueAssignmentScore> scoreUnassignedPos(LMB objFcn, 
			MSSearchProblem[][] objFcnSearch, 
			ArrayList<Integer> unassignedPos) {

		int numSubStates = objFcnSearch[0].length;
		MSSearchProblem complex = objFcnSearch[0][numSubStates-1];

		ArrayList<ResidueAssignment> assignments = new ArrayList<>();

		//get assignments
		for(int splitPos : unassignedPos) {

			if(complex.getNumAssignedPos()==0) {//root, split all possible from unbound
				assignments = getUnboundResidueAssignments(objFcnSearch[0]);
				break;
			}

			else {//can directly score bound state
				assignments.add(getBoundResidueAssignments(objFcnSearch[0], splitPos));
			}
		}

		//assignments.trimToSize();

		//now score each assignment
		ArrayList<ResidueAssignmentScore> assignmentScores = new ArrayList<>();
		for(ResidueAssignment assignment : assignments) {
			BigDecimal score;
			if(computeProduct)
				score = getBoundStateDomainProduct(assignment);
			else
				score = getBoundStateDomainSum(assignment);
			assignmentScores.add(new ResidueAssignmentScore(assignment, score));
		}

		assignmentScores.trimToSize();
		return assignmentScores;
	}

	protected ResidueAssignment getBestResidueAssignment(ArrayList<ResidueAssignmentScore> order) {
		//sort assignments by increasing score
		Collections.sort(order, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore a1, ResidueAssignmentScore a2) {
				return a1.score.compareTo(a2.score)<=0 ? -1 : 1;
			}
		});

		ResidueAssignment best = order.get(0).assignment;
		return best;
	}

	public ArrayList<ResidueAssignmentScore> getAllPossibleAssignments(
			LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores, 
			int numMaxMut) {

		//get number of unassigned positions
		int complex = objFcnSearch[0].length-1;
		ArrayList<Integer> unassignedPos = objFcnSearch[0][complex].getPosNums(false);

		if(unassignedPos.size()==0)
			throw new RuntimeException("ERROR: there are no unassigned positions");

		//score unassigned residues by objfcn
		ArrayList<ResidueAssignmentScore> scores = scoreUnassignedPos(objFcn, objFcnSearch, unassignedPos);

		return scores;
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextResidueAssignment(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores,
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> scores = getAllPossibleAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);		

		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(scores);

		//convert to aas that don't violate the allowed number of mutations
		return getBestAAAssignments(objFcnSearch, best, numMaxMut);
	}

}
