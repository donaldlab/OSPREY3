package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class ResidueOrderDynamicScore extends ResidueOrder {

	public static boolean DEBUG = false;

	private enum ScoreType {
		FSCORE,
		HSCORE,
		DISCREPANCY;
	}

	private ArrayList<ArrayList<ArrayList<BigDecimal>>> residueValues;
	private ArrayList<ArrayList<HashMap<Integer, HashSet<String>>>> residueAAs;
	private BoltzmannCalculator boltzmann;
	private ScoreType scoreType;
	protected HashMap<ResidueAssignment, ArrayList<ArrayList<ArrayList<AAAssignment>>>> ra2AAa;

	public ResidueOrderDynamicScore(MSSearchProblem[][] objFcnSearch, String scoreType) {
		super();
		this.boltzmann = new BoltzmannCalculator();
		this.residueValues = allocateResidueValues(objFcnSearch);
		this.ra2AAa = new HashMap<>();
		this.residueAAs = allocateResidueAAs(objFcnSearch);

		switch(scoreType.toLowerCase()) {
		case "fscore":
			this.scoreType = ScoreType.FSCORE;
			break;
		case "hscore":
			this.scoreType = ScoreType.HSCORE;
			break;
		case "discrepancy":
			this.scoreType = ScoreType.DISCREPANCY;
			break;
		default:
			throw new UnsupportedOperationException("ERROR: unsuported score type: "+scoreType);
		}
	}

	private ArrayList<ArrayList<HashMap<Integer, HashSet<String>>>> allocateResidueAAs(MSSearchProblem[][] objFcnSearch) {
		ArrayList<ArrayList<HashMap<Integer, HashSet<String>>>> ans = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			ans.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				ans.get(state).add(new HashMap<Integer, HashSet<String>>());
				for(int residuePos=0;residuePos<objFcnSearch[state][subState].settings.mutRes.size();++residuePos) {
					ans.get(state).get(subState).put(residuePos, new HashSet<String>());
				}
			}
			ans.get(state).trimToSize();
		}
		ans.trimToSize();
		return ans;
	}

	private void clearResidueAAs() {
		for(int state=0;state<residueAAs.size();++state) {
			for(int subState=0;subState<residueAAs.get(state).size();++subState) {
				for(Integer residuePos : residueAAs.get(state).get(subState).keySet()) {
					residueAAs.get(state).get(subState).get(residuePos).clear();
				}
			}
		}
	}

	private ArrayList<ArrayList<ArrayList<BigDecimal>>> allocateResidueValues(MSSearchProblem[][] objFcnSearch) {
		ArrayList<ArrayList<ArrayList<BigDecimal>>> ans = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			ans.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				ans.get(state).add(new ArrayList<>());
				for(int residuePos=0;residuePos<objFcnSearch[state][subState].settings.mutRes.size();++residuePos) {
					ans.get(state).get(subState).add(BigDecimal.ZERO);
				}
				ans.get(state).get(subState).trimToSize();
			}
			ans.get(state).trimToSize();
		}
		ans.trimToSize();
		return ans;
	}

	protected void setResidueValues(LMB objFcn, MSSearchProblem[][] objFcnSearch, boolean assigned, boolean parallel) {
		ArrayList<ResisueOrderWorker> workers = new ArrayList<>();
		//set all workers
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				boolean isUnbound = subState != objFcnSearch[state].length-1;
				workers.add(new ResisueOrderWorker(objFcnSearch[state][subState], state, subState, isUnbound));
			}
		}

		if(!parallel) {
			for(ResisueOrderWorker w : workers) setResidueValues(objFcn.getCoeffs()[w.state], w.search, w.state, 
					w.subState, w.isUnbound, assigned);
		}
		//execute in parallel
		else {
			workers.parallelStream().forEach(w -> setResidueValues(objFcn.getCoeffs()[w.state], w.search, w.state, 
					w.subState, w.isUnbound, assigned));
		}
	}

	protected boolean isUpperBoundPFProxy(BigDecimal objFcnCoeff, boolean isUnbound) {
		if(objFcnCoeff.compareTo(BigDecimal.ZERO)<0 && !isUnbound ||
				objFcnCoeff.compareTo(BigDecimal.ZERO)>0 && isUnbound)
			return true;
		return false;
	}

	protected void setResidueValues(BigDecimal objFcnCoeff, MSSearchProblem search, int state, int subState, 
			boolean isUnbound, boolean assigned) {

		final boolean useGMECProxy = false;
		double maxVal = search.settings.stericThreshold;
		boolean isUpperBoundPFProxy = isUpperBoundPFProxy(objFcnCoeff, isUnbound);

		ArrayList<Integer> pos1s = search.getNumAssignedPos() > 0 ? search.getPosNums(assigned) : search.getPosNums();
		for(int pos1 : pos1s) {
			double pos1Value = 0;
			//also used pruned rcs at pos?
			ArrayList<Integer> pos1RCs = search.pruneMat.unprunedRCsAtPos(pos1);

			ArrayList<Integer> pos2s = assigned ? pos1s : search.getPosNums();
			//ArrayList<Integer> pos2s = search.getPosNums();
			for(int pos2 : pos2s) {
				if(pos1==pos2) continue;

				if(pos1Value == Double.POSITIVE_INFINITY) continue;

				ArrayList<Integer> pos2RCs = search.pruneMat.unprunedRCsAtPos(pos2);

				// first, find the min pairwise energy over all rc pairs
				double minPairwise = Double.POSITIVE_INFINITY;
				double maxPairwise = Double.NEGATIVE_INFINITY;

				for(int rc1 : pos1RCs) {

					String rc1AAType = search.confSpace.posFlex.get(pos1).RCs.get(rc1).AAType;
					if(!residueAAs.get(state).get(subState).get(pos1).contains(rc1AAType)) continue;

					for(int rc2 : pos2RCs) {

						String rc2AAType = search.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
						if(!residueAAs.get(state).get(subState).get(pos2).contains(rc2AAType)) continue;

						//cap energy at steric threshold
						double pairwise = Math.min(search.emat.getPairwise(pos1, rc1, pos2, rc2), maxVal);

						minPairwise = Math.min(minPairwise, pairwise);
						maxPairwise = Math.max(maxPairwise, pairwise);

					}
				}


				//adjust minpairwise upwards for lowerbound case
				if(!isUpperBoundPFProxy) {
					if(minPairwise != maxVal && maxPairwise == maxVal)
						minPairwise = (minPairwise + maxPairwise)/2.0;
						maxPairwise = minPairwise;
				}

				if(useGMECProxy) {
					//TRY MIN/MAX ENERGIES
					if(isUpperBoundPFProxy) 
						pos1Value += minPairwise;
					else
						pos1Value += maxPairwise;
				} else {
					//TRY HARMONIC MEANS
					//computing lower bound partition function proxy
					double pos2Value = 0;
					for(int rc1 : pos1RCs) {

						String rc1AAType = search.confSpace.posFlex.get(pos1).RCs.get(rc1).AAType;
						if(!residueAAs.get(state).get(subState).get(pos1).contains(rc1AAType)) continue;

						for(int rc2 : pos2RCs) {

							String rc2AAType = search.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
							if(!residueAAs.get(state).get(subState).get(pos2).contains(rc2AAType)) continue;

							//cap energy at steric threshold
							double pairwise = Math.min(search.emat.getPairwise(pos1, rc1, pos2, rc2), maxVal);
							double normalizedPairwise = 0;

							if(isUpperBoundPFProxy) {
								normalizedPairwise = pairwise - minPairwise;
							}

							else {
								//adjust pairwise so that it's >= min pairwise
								pairwise = Math.max(pairwise, minPairwise);
								normalizedPairwise = pairwise - minPairwise;
							}

							if (normalizedPairwise != 0) {
								pos2Value += 1.0/normalizedPairwise;
							}
						}
					}

					int pos1Size = search.countRcsAtPosForAAs(search.pruneMat, pos1, residueAAs.get(state).get(subState).get(pos1), false);
					int pos2Size = search.countRcsAtPosForAAs(search.pruneMat, pos2, residueAAs.get(state).get(subState).get(pos2), false);

					pos2Value = (pos1Size*pos2Size - 1)/pos2Value;
					if(Double.isNaN(pos2Value) || Double.isInfinite(pos2Value)) continue;
					pos1Value += pos2Value;
				}

			}

			if(Double.isNaN(pos1Value)) pos1Value = 0;
			BigDecimal val = boltzmann.calc(pos1Value).setScale(64, RoundingMode.HALF_UP);

			//val can be 0 for bound state but not for unbound states
			if(isUnbound) {
				val = val.compareTo(BigDecimal.ZERO)==0 ? BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP) : val;
			}

			/*
			//multiply by number of rotamers for the upper bound partition function proxy
			if(isUpperBoundPFProxy) {
				//multiply by the number of rotamers at pos
				//val = val.multiply(BigDecimal.valueOf(pos1Size));
				int numAAsAtPos = residueAAs.get(state).get(subState).get(pos1).size();
				int pos1Size = search.countRcsAtPosForAAs(search.pruneMat, pos1, residueAAs.get(state).get(subState).get(pos1), false);
				val = val.multiply(BigDecimal.valueOf((double)pos1Size/(double)numAAsAtPos));
			}
			*/

			residueValues.get(state).get(subState).set(pos1, val);
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

	protected ArrayList<ResidueAssignment> getUnboundResidueAssignments(MSSearchProblem[] search) {
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

	protected ResidueAssignment getBoundResidueAssignments(MSSearchProblem[] search, 
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

	private BigDecimal getStateScoreHScore(int state, ResidueAssignment assignment) {		
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		//unbound states
		for(int subState=0;subState<assignment.length()-1;++subState) {
			ArrayList<Integer> unboundPos = assignment.get(subState);
			if(unboundPos.size()==0) continue;
			if(unboundPos.size()>1) throw new RuntimeException("ERROR: unbound state was split into more than one position");
			int pos = unboundPos.get(0);//contains at most one value		
			ans = ans.multiply(residueValues.get(state).get(subState).get(pos));
		}
		if(ans.compareTo(BigDecimal.ZERO)==0) ans = new BigDecimal("1e-64").setScale(64, RoundingMode.HALF_UP);

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

	private BigDecimal getStateScoreFScore1(int state, 
			ResidueAssignment assignment,
			KStarScore gScore) {
		ArrayList<BigDecimal> stateFScores = new ArrayList<>();

		//iterate through substates
		for(int subState=0;subState<assignment.length();++subState) {
			//get g score from previous assignment (will be 0 if root)
			BigDecimal g = gScore.getPartitionFunction(subState) == null ? 
					BigDecimal.ZERO : gScore.getPartitionFunction(subState).getValues().qstar;
			g = g.setScale(64, RoundingMode.HALF_UP);

			//add h score to gscore
			ArrayList<Integer> unboundPos = assignment.get(subState);
			BigDecimal h = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			for(int pos : unboundPos) h = h.add(residueValues.get(state).get(subState).get(pos));

			//add up g and h scores for each substate
			stateFScores.add(g.add(h));
		}

		//then do final division
		BigDecimal denom = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int subState=0;subState<stateFScores.size()-1;++subState) denom = denom.multiply(stateFScores.get(subState));
		if(denom.compareTo(BigDecimal.ZERO)==0) denom = new BigDecimal("1e-64").setScale(64, RoundingMode.HALF_UP);

		int complex = stateFScores.size()-1;
		BigDecimal numer = stateFScores.get(complex).setScale(64, RoundingMode.HALF_UP);

		BigDecimal ans = numer.divide(denom, RoundingMode.HALF_UP);
		return ans;
	}

	private BigDecimal getStateScoreFScorePartial(int state, 
			MSSearchProblem[][] objFcnSearch, 
			ResidueAssignment assignment) {

		ArrayList<BigDecimal> stateFScores = new ArrayList<>();

		//iterate through substates
		for(int subState=0;subState<assignment.length();++subState) {
			//get g score from previous assignment (will be 0 if root)
			BigDecimal g = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ArrayList<Integer> assignedPos = objFcnSearch[state][subState].getPosNums(true); 
			for(int pos : assignedPos) g = g.add(residueValues.get(state).get(subState).get(pos));

			//add h score to gscore
			ArrayList<Integer> unboundPos = assignment.get(subState);
			BigDecimal h = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			for(int pos : unboundPos) h = h.add(residueValues.get(state).get(subState).get(pos));

			//add up g and h scores for each substate
			stateFScores.add(g.add(h));
		}

		//then do final division
		BigDecimal denom = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int subState=0;subState<stateFScores.size()-1;++subState) denom = denom.multiply(stateFScores.get(subState));
		if(denom.compareTo(BigDecimal.ZERO)==0) denom = new BigDecimal("1e-64").setScale(64, RoundingMode.HALF_UP);

		int complex = stateFScores.size()-1;
		BigDecimal numer = stateFScores.get(complex).setScale(64, RoundingMode.HALF_UP);

		BigDecimal ans = numer.divide(denom, RoundingMode.HALF_UP);
		return ans;
	}

	private BigDecimal getStateScoreGScore(int state, 
			MSSearchProblem[][] objFcnSearch, 
			ResidueAssignment assignment) {

		ArrayList<BigDecimal> stateFScores = new ArrayList<>();

		//iterate through substates
		for(int subState=0;subState<assignment.length();++subState) {
			//get g score from previous assignment (will be 0 if root)
			BigDecimal g = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ArrayList<Integer> assignedPos = objFcnSearch[state][subState].getPosNums(true);
			for(int pos : assignedPos) g = g.add(residueValues.get(state).get(subState).get(pos));

			//add up g and h scores for each substate
			stateFScores.add(g);
		}

		//then do final division
		BigDecimal denom = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int subState=0;subState<stateFScores.size()-1;++subState) denom = denom.multiply(stateFScores.get(subState));
		if(denom.compareTo(BigDecimal.ZERO)==0) denom = new BigDecimal("1e-64").setScale(64, RoundingMode.HALF_UP);

		int complex = stateFScores.size()-1;
		BigDecimal numer = stateFScores.get(complex).setScale(64, RoundingMode.HALF_UP);

		BigDecimal ans = numer.divide(denom, RoundingMode.HALF_UP);
		return ans;
	}

	private BigDecimal getStateScoreFScoreFull(int state, 
			MSSearchProblem[][] objFcnSearch, 
			ResidueAssignment assignment) {

		ArrayList<BigDecimal> stateFScores = new ArrayList<>();

		//iterate through substates
		for(int subState=0;subState<assignment.length();++subState) {
			//get g score from previous assignment (will be 0 if root)
			BigDecimal g = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ArrayList<Integer> assignedPos = objFcnSearch[state][subState].getPosNums(true);
			for(int pos : assignedPos) g = g.add(residueValues.get(state).get(subState).get(pos));

			//add h score to gscore
			BigDecimal h = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ArrayList<Integer> unassignedPos = objFcnSearch[state][subState].getPosNums(false);
			for(int pos : unassignedPos) h = h.add(residueValues.get(state).get(subState).get(pos));

			//add up g and h scores for each substate
			stateFScores.add(g.add(h));
		}

		//then do final division
		BigDecimal denom = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		for(int subState=0;subState<stateFScores.size()-1;++subState) denom = denom.multiply(stateFScores.get(subState));
		if(denom.compareTo(BigDecimal.ZERO)==0) denom = new BigDecimal("1e-64").setScale(64, RoundingMode.HALF_UP);

		int complex = stateFScores.size()-1;
		BigDecimal numer = stateFScores.get(complex).setScale(64, RoundingMode.HALF_UP);

		BigDecimal ans = numer.divide(denom, RoundingMode.HALF_UP);
		return ans;
	}

	//includes k* scores from current node
	private BigDecimal getScore(ResidueAssignment assignment, 
			MSSearchProblem[][] objFcnSearch,
			BigDecimal[] coeffs,
			BigDecimal fScore) {

		int numStates = coeffs.length;
		BigDecimal ans = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);

		//get assignment score
		for(int state=0;state<numStates;++state) {
			//sign of 0 does not contribute to score
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;

			if(scoreType == ScoreType.DISCREPANCY) {
				BigDecimal fScoreFull = getStateScoreFScoreFull(state, objFcnSearch, assignment);
				BigDecimal fScorePartial = getStateScoreFScorePartial(state, objFcnSearch, assignment);
				ans = ans.add(fScoreFull.subtract(fScorePartial).abs().multiply(BigDecimal.valueOf(-1)));
			}
			else if(scoreType == ScoreType.FSCORE)
				ans = ans.add(coeffs[state].multiply(getStateScoreFScorePartial(state, objFcnSearch, assignment)));

			else if(scoreType == ScoreType.HSCORE)
				ans = ans.add(coeffs[state].multiply(getStateScoreHScore(state, assignment)));

			else
				throw new UnsupportedOperationException("ERROR: scoreType must be one of {DISCREPANCY, FSCORE, HSCORE}");
		}

		return ans;
	}

	protected ArrayList<ResidueAssignmentScore> scoreUnassignedPos(LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores) {

		BigDecimal[] stateVals = new BigDecimal[objFcnScores.length];
		for(int i=0;i<stateVals.length;++i) stateVals[i] = objFcnScores[i].getScore();
		BigDecimal fScore = objFcn.eval(stateVals);

		//now score each assignment
		ArrayList<ResidueAssignmentScore> assignmentScores = new ArrayList<>();
		for(ResidueAssignment residueAssignment : ra2AAa.keySet()) {
			//if coeff[state]<0: want ub ratio, so smallest unbound state, largest bound state
			//if coeff[state]>0: want lb ratio, so largest unbound state, smallest bound state
			BigDecimal score = getScore(residueAssignment, objFcnSearch, objFcn.getCoeffs(), fScore);
			assignmentScores.add(new ResidueAssignmentScore(residueAssignment, score));
		}

		//assignmentScores.trimToSize();
		return assignmentScores;
	}

	protected ResidueAssignment getBestResidueAssignment(ArrayList<ResidueAssignmentScore> order) {
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

	protected ArrayList<ArrayList<ArrayList<AAAssignment>>> getAllowedAAAsignments(
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

	public ArrayList<ResidueAssignmentScore> scoreAllResidueAssignments(
			LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores, 
			int numMaxMut
			) {
		//get number of unassigned positions
		int complex = objFcnSearch[0].length-1;
		ArrayList<Integer> unassignedPos = objFcnSearch[0][complex].getPosNums(false);

		if(unassignedPos.size()==0)
			throw new RuntimeException("ERROR: there are no unassigned positions");

		clearResidue2AAAssignments();
		clearResidueAAs();
		computeAllowedAAsAtPos(objFcnSearch, unassignedPos, numMaxMut);//add unassigned

		//"g-score": value assigned residues
		setResidueValues(objFcn, objFcnSearch, true, MSKStarNode.PARALLEL_EXPANSION);

		//"h-score": value unassigned residues
		setResidueValues(objFcn, objFcnSearch, false, MSKStarNode.PARALLEL_EXPANSION);

		//score unassigned residues by objfcn
		ArrayList<ResidueAssignmentScore> scores = scoreUnassignedPos(objFcn, objFcnSearch, objFcnScores);

		return scores;
	}

	private void computeAllowedAAsAtPos(MSSearchProblem[][] objFcnSearch, ArrayList<Integer> unassignedPos, int numMaxMut) {
		//compute for assigned pos
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int substate=0;substate<objFcnSearch[state].length;++substate) {
				for(int pos : objFcnSearch[state][substate].getPosNums(true)) {

					MSSearchProblem search = objFcnSearch[state][substate];

					if(search.settings.AATypeOptions.get(pos).size() != 1)
						throw new RuntimeException("ERROR: there must be exactly one AA per residue assignment");

					String allowedAA = search.settings.AATypeOptions.get(pos).get(0);

					residueAAs.get(state).get(substate).get(pos).add(allowedAA);
				}
			}
		}

		//compute for unassigned pos
		int numSubStates = objFcnSearch[0].length;
		MSSearchProblem complex = objFcnSearch[0][numSubStates-1];
		ArrayList<ResidueAssignment> residueAssignments = new ArrayList<>();
		//get assignments
		for(int splitPos : unassignedPos) {

			if(complex.getNumAssignedPos()==0) {//root, split >=2 pos
				residueAssignments = getUnboundResidueAssignments(objFcnSearch[0]);
				break;
			}

			else {//can directly score bound state
				residueAssignments.add(getBoundResidueAssignments(objFcnSearch[0], splitPos));
			}
		}

		//now, set allowed AAs from residue assignments
		for(ResidueAssignment residueAssignment : residueAssignments) {
			ArrayList<ArrayList<ArrayList<AAAssignment>>> aaAssignments = getAllowedAAAsignments(objFcnSearch, residueAssignment, numMaxMut);
			storeResidue2AAAssignments(residueAssignment, aaAssignments);

			int numSplits = aaAssignments.get(0).size();

			for(int state=0;state<objFcnSearch.length;++state) {
				for(int split=0;split<numSplits;++split) {
					for(int substate=0;substate<residueAssignment.length();++substate) {
						ArrayList<AAAssignment> as = aaAssignments.get(substate).get(split);
						for(AAAssignment a : as) {
							int rPos = a.residuePos;
							int aPos = a.AATypePos;
							HashSet<String> aasAtPos = residueAAs.get(state).get(substate).get(rPos);
							String allowedAA = objFcnSearch[state][substate].settings.AATypeOptions.get(rPos).get(aPos);
							if(!aasAtPos.contains(allowedAA))
								aasAtPos.add(allowedAA);
						}
					}
				}
			}
		}

		if(DEBUG) {
			//do sanity check to make sure that
			//1) assigned residue only contains 1 single AA
			for(int state=0;state<objFcnSearch.length;++state) {
				for(int substate=0;substate<objFcnSearch[state].length;++substate) {
					for(int pos : objFcnSearch[state][substate].getPosNums(true)) {
						HashSet<String> aasAtPos = residueAAs.get(state).get(substate).get(pos);
						if(aasAtPos.size() != 1)
							throw new RuntimeException("ERROR: there must be exactly one AA per residue assignment");
					}
				}
			}

			//2) each aa is contained in AA options in the search.settings
			for(int state=0;state<objFcnSearch.length;++state) {
				for(int substate=0;substate<objFcnSearch[state].length;++substate) {
					for(int pos : objFcnSearch[state][substate].getPosNums(false)) {
						HashSet<String> aasAtPos = residueAAs.get(state).get(substate).get(pos);
						ArrayList<String> allowedAAsAtPos = objFcnSearch[state][substate].settings.AATypeOptions.get(pos);
						for(String aa : aasAtPos) {
							if(!allowedAAsAtPos.contains(aa))
								throw new RuntimeException("ERROR: "+aa+" is not allowed");
						}
					}
				}
			}
		}
	}

	protected void storeResidue2AAAssignments(ResidueAssignment residueAssignment, 
			ArrayList<ArrayList<ArrayList<AAAssignment>>> aaAssignments) {
		ra2AAa.put(residueAssignment, aaAssignments);
	}

	protected ArrayList<ArrayList<ArrayList<AAAssignment>>> getAAAssignments(ResidueAssignment residueAssignment) {
		ArrayList<ArrayList<ArrayList<AAAssignment>>> ans = ra2AAa.get(residueAssignment);
		return ans;
	}

	protected void clearResidue2AAAssignments() {
		ra2AAa.clear();
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores,
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> scores = scoreAllResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);

		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(scores);

		//convert to aas that don't violate the allowed number of mutations
		return getAAAssignments(best);
	}

}
