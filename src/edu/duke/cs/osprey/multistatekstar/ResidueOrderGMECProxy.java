package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings("serial")
public class ResidueOrderGMECProxy extends ResidueOrderGMEC {

	private class ResisueOrderWorker {
		public int state;
		public int subState;
		public MSSearchProblem search;
		boolean isUnbound;
		
		public ResisueOrderWorker(MSSearchProblem search, int state, int subState, boolean isUnbound) {
			this.search = search;
			this.state = state;
			this.subState = subState;
			this.isUnbound = isUnbound;
		}
	}
	
	private enum ScoreType {
		HMEAN,
		LOWEST;
	}
	
	public static boolean DEBUG = false;
	private ScoreType scoreType;

	public ResidueOrderGMECProxy(MSSearchProblem[][] objFcnSearch, String method) {
		super(objFcnSearch);
		this.scoreType = getScoreType(method);
	}
	
	private ScoreType getScoreType(String method) {
		switch(method.toLowerCase()) {
		case "hmean":
			return ScoreType.HMEAN;
		case "lowest":
			return ScoreType.LOWEST;
		default:
			throw new UnsupportedOperationException("ERROR: method must be in the set {hmean, lowest}");
		}
	}

	protected void computeStateScores(LMB objFcn, MSSearchProblem[][] objFcnSearch, boolean assigned, boolean parallel) {
		ArrayList<ResisueOrderWorker> workers = new ArrayList<>();
		//set all workers
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				boolean isUnbound = subState != objFcnSearch[state].length-1;
				workers.add(new ResisueOrderWorker(objFcnSearch[state][subState], state, subState, isUnbound));
			}
		}

		if(!parallel) {
			for(ResisueOrderWorker w : workers) computeStateScores(objFcn.getCoeffs()[w.state], w.search, w.state, 
					w.subState, w.isUnbound, assigned);
		}
		//execute in parallel
		else {
			workers.parallelStream().forEach(w -> computeStateScores(objFcn.getCoeffs()[w.state], w.search, w.state, 
					w.subState, w.isUnbound, assigned));
		}
	}

	protected boolean isUpperBoundPFProxy(BigDecimal objFcnCoeff, boolean isUnbound) {
		if(objFcnCoeff.compareTo(BigDecimal.ZERO)<0 && !isUnbound ||
				objFcnCoeff.compareTo(BigDecimal.ZERO)>0 && isUnbound)
			return true;
		return false;
	}

	protected void computeStateScores(BigDecimal objFcnCoeff, MSSearchProblem search, int state, int subState, 
			boolean isUnbound, boolean assigned) {

		double maxVal = search.settings.stericThreshold;
		//boolean isUpperBoundPFProxy = isUpperBoundPFProxy(objFcnCoeff, isUnbound);

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
					if(!pos2AAs.get(state).get(subState).get(pos1).contains(rc1AAType)) continue;

					for(int rc2 : pos2RCs) {

						String rc2AAType = search.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
						if(!pos2AAs.get(state).get(subState).get(pos2).contains(rc2AAType)) continue;

						//cap energy at steric threshold
						double pairwise = Math.min(search.emat.getPairwise(pos1, rc1, pos2, rc2), maxVal);

						minPairwise = Math.min(minPairwise, pairwise);
						maxPairwise = Math.max(maxPairwise, pairwise);
					}
				}

				//adjust minpairwise upwards for lowerbound case
				/*
				if(isUpperBoundPFProxy) {
					if(minPairwise != maxVal && maxPairwise == maxVal) {
						maxPairwise = (minPairwise + maxPairwise)/2.0;
						minPairwise = maxPairwise;
					}
				}
				*/
				
				if(scoreType == ScoreType.LOWEST) {
					/*
					//TRY MIN/MAX ENERGIES
					if(isUpperBoundPFProxy)
						pos1Value += maxPairwise;
					else
						pos1Value += minPairwise;
					*/
					pos1Value += minPairwise;
				} else if(scoreType == ScoreType.HMEAN) {
					//TRY HARMONIC MEANS
					//computing lower bound partition function proxy
					double pos2Value = 0;
					for(int rc1 : pos1RCs) {

						String rc1AAType = search.confSpace.posFlex.get(pos1).RCs.get(rc1).AAType;
						if(!pos2AAs.get(state).get(subState).get(pos1).contains(rc1AAType)) continue;

						for(int rc2 : pos2RCs) {

							String rc2AAType = search.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
							if(!pos2AAs.get(state).get(subState).get(pos2).contains(rc2AAType)) continue;

							//cap energy at steric threshold
							double pairwise = Math.min(search.emat.getPairwise(pos1, rc1, pos2, rc2), maxVal);
							//if(isUpperBoundPFProxy) pairwise = Math.max(pairwise, minPairwise);
							
							double normalizedPairwise = pairwise - minPairwise;

							if (normalizedPairwise != 0) {
								pos2Value += 1.0/normalizedPairwise;
							}
						}
					}

					int pos1Size = search.countRcsAtPosForAAs(search.pruneMat, pos1, pos2AAs.get(state).get(subState).get(pos1), false);
					int pos2Size = search.countRcsAtPosForAAs(search.pruneMat, pos2, pos2AAs.get(state).get(subState).get(pos2), false);

					pos2Value = (pos1Size*pos2Size - 1)/pos2Value;
					if(Double.isNaN(pos2Value) || Double.isInfinite(pos2Value)) continue;
					pos1Value += pos2Value;
				}

			}

			if(Double.isNaN(pos1Value)) pos1Value = 0;
			pos1Value = Math.min(pos1Value, maxVal);

			pos2Score.get(state).get(subState).set(pos1, pos1Value);
		}
	}

	protected BigDecimal getStateScoreGScore(int state, 
			MSSearchProblem[][] objFcnSearch, 
			ResidueAssignment assignment) {

		ArrayList<BigDecimal> stateFScores = new ArrayList<>();

		//iterate through substates
		for(int subState=0;subState<assignment.length();++subState) {
			//get g score from previous assignment (will be 0 if root)
			BigDecimal g = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
			ArrayList<Integer> assignedPos = objFcnSearch[state][subState].getPosNums(true);
			for(int pos : assignedPos) g = g.add(BigDecimal.valueOf(pos2Score.get(state).get(subState).get(pos)));

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

	protected ResidueAssignment getBestResidueAssignment(ArrayList<ResidueAssignmentScore> order) {
		//sort assignments by decreasing score
		Collections.sort(order, new Comparator<ResidueAssignmentScore>() {
			@Override
			public int compare(ResidueAssignmentScore a1, ResidueAssignmentScore a2) {
				return a1.score.compareTo(a2.score)<=0 ? -1 : 1;
			}
		});

		ResidueAssignment best = order.get(0).assignment;
		return best;
	}
	
	public ArrayList<ResidueAssignmentScore> scoreResidueAssignments(
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

		clearStoredStructures();
		
		computeAAsAndResAssignments(objFcnSearch, unassignedPos, numMaxMut);

		//"g-score": value assigned residues
		computeStateScores(objFcn, objFcnSearch, true, MSKStarNode.PARALLEL_EXPANSION);

		//"h-score": value unassigned residues
		computeStateScores(objFcn, objFcnSearch, false, MSKStarNode.PARALLEL_EXPANSION);
		
		//score unassigned residues by objfcn
		ArrayList<ResidueAssignmentScore> scores = scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores);

		return scores;
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(
			LMB objFcn,
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores,
			int numMaxMut
			) {
		ArrayList<ResidueAssignmentScore> scores = scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);

		//order unassigned residues by score and return best residue
		ResidueAssignment best = getBestResidueAssignment(scores);

		//convert to aas that don't violate the allowed number of mutations
		return getAAAssignments(best);
	}

}
