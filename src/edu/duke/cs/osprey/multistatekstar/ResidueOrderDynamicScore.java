package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class ResidueOrderDynamicScore extends ResidueOrder {

	private ArrayList<ArrayList<ArrayList<BigDecimal>>> residueValues;
	private BoltzmannCalculator boltzmann;

	public ResidueOrderDynamicScore(MSSearchProblem[][] objFcnSearch) {
		super();
		this.boltzmann = new BoltzmannCalculator();
		this.residueValues = allocate(objFcnSearch);
	}

	private ArrayList<ArrayList<ArrayList<BigDecimal>>> allocate(MSSearchProblem[][] objFcnSearch) {
		//create a priority queue for each state and substate
		ArrayList<ArrayList<ArrayList<BigDecimal>>> residueValues = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			residueValues.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				residueValues.get(state).add(new ArrayList<>());
				for(int residuePos=0;residuePos<objFcnSearch[state][subState].settings.mutRes.size();++residuePos) {
					residueValues.get(state).get(subState).add(BigDecimal.ZERO);
				}
				residueValues.get(state).get(subState).trimToSize();
			}
			residueValues.get(state).trimToSize();
		}
		residueValues.trimToSize();
		return residueValues;
	}

	protected void setResidueValues(MSSearchProblem[][] objFcnSearch, boolean assigned, boolean parallel) {
		ArrayList<ResisueOrderWorker> workers = new ArrayList<>();
		//set all workers
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				workers.add(new ResisueOrderWorker(objFcnSearch[state][subState], state, subState));
			}
		}

		if(!parallel) {
			for(ResisueOrderWorker w : workers) setResidueValues(w.search, w.state, w.subState, w.search.getPosNums(assigned));
		}
		//execute in parallel
		else {
			workers.parallelStream().forEach(w -> setResidueValues(w.search, w.state, w.subState, w.search.getPosNums(assigned)));
		}
	}

	protected void setResidueValues(MSSearchProblem search, int state, int subState, ArrayList<Integer> positions) {

		int numPos = search.getNumPos();

		for(int pos1 : positions) {
			double pos1Value = 0;
			//also used pruned rcs at pos?
			ArrayList<Integer> pos1RCs = search.pruneMat.unprunedRCsAtPos(pos1);

			for(int pos2=0;pos2<numPos;++pos2) {
				if(pos1==pos2) continue;

				ArrayList<Integer> pos2RCs = search.pruneMat.unprunedRCsAtPos(pos2);

				// first, find the min pairwise energy over all rc pairs
				double minPairwise = Double.POSITIVE_INFINITY;				
				for(int rc1 : pos1RCs) {
					for(int rc2 : pos2RCs) {
						minPairwise = Math.min(minPairwise, search.emat.getPairwise(pos1, rc1, pos2, rc2));
					}
				}

				//compute the harmonic average of normalized pairwise energies
				double pos2Value = 0;
				for(int rc1 : pos1RCs) {
					for(int rc2 : pos2RCs) {
						double normalizedPairwise = search.emat.getPairwise(pos1, rc1, pos2, rc2) - minPairwise;
						if (normalizedPairwise != 0) {
							pos2Value += 1.0/normalizedPairwise;
						}
					}
				}

				pos2Value = (pos1RCs.size()*pos2RCs.size() - 1)/pos2Value;
				pos1Value += pos2Value;
			}

			if(Double.isNaN(pos1Value)) pos1Value = 0;
			residueValues.get(state).get(subState).set(pos1, (boltzmann.calc(pos1Value)).setScale(64, RoundingMode.HALF_UP));
		}
	}

	private void generatePermutations(MSSearchProblem[] search,
			ArrayList<ArrayList<String>> input, ArrayList<ResidueAssignment> output, 
			ArrayList<String> buf, int depth) {

		if(depth==input.size()) {//each unbound state is assigned
			ArrayList<ArrayList<Integer>> assignments = new ArrayList<>();
			ArrayList<Integer> complex = new ArrayList<>();

			for(int subState=0;subState<depth;++subState) {
				assignments.add(new ArrayList<>());
				assignments.get(subState).add(search[subState].flexRes.indexOf(buf.get(subState)));

				complex.add(search[search.length-1].flexRes.indexOf(buf.get(subState)));
				assignments.get(subState).trimToSize();
			}

			complex.trimToSize();
			assignments.add(complex);
			assignments.trimToSize();

			output.add(new ResidueAssignment(assignments));
			return;
		}

		for(int i=0;i<input.get(depth).size();++i) {
			buf.set(depth, input.get(depth).get(i));
			generatePermutations(search, input, output, buf, depth+1);
		}
	}

	private ArrayList<ResidueAssignment> getUnboundResidueAssignments(MSSearchProblem[] search) {
		ArrayList<ArrayList<String>> unassignedRes = new ArrayList<>();
		ArrayList<String> buf = new ArrayList<>();
		//get unbound states unassigned flex res
		for(int subState=0;subState<search.length-1;++subState) {
			unassignedRes.add(search[subState].getResidues(false));
			buf.add("");
		}
		unassignedRes.trimToSize();
		buf.trimToSize();

		ArrayList<ResidueAssignment> ans = new ArrayList<>();
		//all n combinations of unbound state, taking 1 residue from each unbound state
		generatePermutations(search, unassignedRes, ans, buf, 0);		
		ans.trimToSize();
		return ans;
	}

	private ResidueAssignment getBoundResidueAssignments(MSSearchProblem[] search, 
			int splitPos) {

		ArrayList<Integer> complex = new ArrayList<>(); 
		complex.add(splitPos); complex.trimToSize();

		String complexRes = search[search.length-1].flexRes.get(splitPos);

		ArrayList<ArrayList<Integer>> assignments = new ArrayList<>();
		for(int subState=0;subState<search.length-1;++subState) {
			assignments.add(new ArrayList<>());
			int pos = search[subState].flexRes.indexOf(complexRes);
			if(pos != -1) 
				assignments.get(subState).add(pos);
			assignments.get(subState).trimToSize();
		}

		assignments.add(complex); assignments.trimToSize();

		ResidueAssignment ans = new ResidueAssignment(assignments);
		return ans;
	}

	private BigDecimal getStateScore(int state, ResidueAssignment assignment) {
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		//unbound states
		for(int subState=0;subState<assignment.length()-1;++subState) {
			ArrayList<Integer> unboundPos = assignment.get(subState);
			if(unboundPos.size()==0) continue;
			if(unboundPos.size()>1) throw new RuntimeException("ERROR: unbound state was split into more than one position");
			int pos = unboundPos.get(0);//contains at most one value
			ans = ans.divide(residueValues.get(state).get(subState).get(pos), RoundingMode.HALF_UP);
		}

		//bound state
		int subState = assignment.length()-1;
		ArrayList<Integer> boundPos = assignment.get(subState);
		BigDecimal numerator = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		for(int pos : boundPos) {
			numerator = numerator.add(residueValues.get(state).get(subState).get(pos));
		}

		ans = numerator.divide(ans, RoundingMode.HALF_UP);
		return ans;
	}

	private BigDecimal getFScore(ResidueAssignment assignment, BigDecimal[] coeffs) {
		int numStates = coeffs.length;
		BigDecimal ans = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		// get assignment score
		for(int state=0;state<numStates;++state) {
			//sign of 0 does not contribute to score
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			ans = ans.add(coeffs[state].multiply(getStateScore(state, assignment)));
		}

		return ans;
	}

	private ArrayList<ResidueAssignmentScore> scoreUnassignedPos(LMB objFcn, MSSearchProblem[][] objFcnSearch, 
			ArrayList<Integer> unassignedPos) {

		int numSubStates = objFcnSearch[0].length;
		MSSearchProblem complex = objFcnSearch[0][numSubStates-1];

		ArrayList<ResidueAssignment> assignments = new ArrayList<>();

		//get assignments
		for(int splitPos : unassignedPos) {

			if(complex.getNumAssignedPos()==0) {//root
				assignments = getUnboundResidueAssignments(objFcnSearch[0]);
				break;
			}

			else {//can directly score bound state
				assignments.add(getBoundResidueAssignments(objFcnSearch[0], splitPos));
			}
		}

		assignments.trimToSize();

		//if coeff[state]<0: want ub ratio, so smallest unbound state, largest bound state
		//if coeff[state]>0: want lb ratio, so largest unbound state, smallest bound state
		BigDecimal[] coeffs = objFcn.getCoeffs();
		ArrayList<ResidueAssignmentScore> assignmentScores = new ArrayList<>();
		for(ResidueAssignment assignment : assignments) {
			BigDecimal score = getFScore(assignment, coeffs);
			assignmentScores.add(new ResidueAssignmentScore(assignment, score));
		}

		assignmentScores.trimToSize();
		return assignmentScores;
	}

	private ResidueAssignment getBestResidueAssignment(ArrayList<ResidueAssignmentScore> order) {
		//sort assignments by decreasing score
		Collections.sort(order, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore a1, ResidueAssignmentScore a2) {
				return a1.score.compareTo(a2.score)>=0 ? -1 : 1;
			}
		});

		ResidueAssignment best = order.get(0).assignment;
		return best;
	}

	/**
	 * get bound state aa assignments that don't exceed the number of allowed 
	 * mutations from wt
	 * @param search
	 * @param splitPos
	 * @param numMaxMut
	 * @return
	 */
	private ArrayList<ArrayList<AAAssignment>> getBoundAAAssignments(
			MSSearchProblem[][] objFcnSearch, 
			ArrayList<Integer> splitPos, 
			int numMaxMut
			) {
		MSSearchProblem complex = objFcnSearch[0][objFcnSearch[0].length-1];

		ArrayList<ArrayList<AAAssignment>> ans = new ArrayList<>();
		String[] wt = MSKStarNode.WT_SEQS.get(0);
		String[] buf = new String[wt.length];
		getBoundAAAssignmentsHelper(complex.settings.AATypeOptions, ans, splitPos, wt, buf, 0, 0, numMaxMut);

		ans.trimToSize();
		return ans;
	}

	private ArrayList<ArrayList<ArrayList<AAAssignment>>> getBestAAAssignments(
			MSSearchProblem[][] objFcnSearch,
			ResidueAssignment residueAssignment, 
			int numMaxMut
			) {

		//we only care about bound state assignment
		ArrayList<Integer> splitPos = residueAssignment.get(residueAssignment.length()-1);
		ArrayList<ArrayList<AAAssignment>> boundAAAssignments = getBoundAAAssignments(objFcnSearch, splitPos, numMaxMut);

		ArrayList<ArrayList<ArrayList<AAAssignment>>> ans = new ArrayList<>();
		int state = 0;
		int numSubStates = objFcnSearch[state].length;
		MSSearchProblem boundState = objFcnSearch[state][numSubStates-1];
		//get corresponding unbound state splits
		for(int subState=0;subState<numSubStates-1;++subState) {
			MSSearchProblem unbound = objFcnSearch[state][subState];
			ans.add(getUnboundAAAssignments(boundState, boundAAAssignments, unbound));
		}
		ans.add(boundAAAssignments);

		ans.trimToSize();
		return ans;
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextResidueAssignment(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch, 
			int numMaxMut
			) {
		//get number of unassigned positions
		int complex = objFcnSearch[0].length-1;
		ArrayList<Integer> unassignedPos = objFcnSearch[0][complex].getPosNums(false);

		if(unassignedPos.size()==0)
			throw new RuntimeException("ERROR: there are no unassigned positions");

		if(unassignedPos.size() > 1) {
			//"g-score": value assigned residues
			setResidueValues(objFcnSearch, true, MSKStarNode.PARALLEL_EXPANSION);

			//"h-score": value unassigned residues
			setResidueValues(objFcnSearch, false, MSKStarNode.PARALLEL_EXPANSION);
		}

		//score unassigned residues by objfcn
		ArrayList<ResidueAssignmentScore> scores = scoreUnassignedPos(objFcn, objFcnSearch, unassignedPos);

		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(scores);

		//convert to aas that don't violate the allowed number of mutations
		return getBestAAAssignments(objFcnSearch, best, numMaxMut);
	}
}
