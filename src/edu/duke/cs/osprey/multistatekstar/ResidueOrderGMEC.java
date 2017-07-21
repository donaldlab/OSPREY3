package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfSearchFactory;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.structure.Residue;

@SuppressWarnings("serial")
public class ResidueOrderGMEC extends ResidueOrder {

	private class ResisueOrderWorker {
		public int state;
		public int subState;
		boolean isUnbound;

		public ResisueOrderWorker(int state, int subState, boolean isUnbound) {
			this.state = state;
			this.subState = subState;
			this.isUnbound = isUnbound;
		}
	}

	public static boolean DEBUG = false;

	protected ArrayList<ArrayList<ArrayList<Double>>> pos2Score;
	protected ArrayList<ArrayList<HashMap<Integer, HashSet<String>>>> pos2AAs;
	protected HashMap<ResidueAssignment, ArrayList<ArrayList<ArrayList<AAAssignment>>>> ra2AAa;
	
	protected BoltzmannCalculator boltzmann;

	public ResidueOrderGMEC(MSSearchProblem[][] objFcnSearch) {
		super();
		this.pos2Score = allocatePos2Score(objFcnSearch);
		this.pos2AAs = allocatePos2AAs(objFcnSearch);
		this.ra2AAa = new HashMap<>();
		
		this.boltzmann = new BoltzmannCalculator();
	}

	private ArrayList<ArrayList<ArrayList<Double>>> allocatePos2Score(MSSearchProblem[][] objFcnSearch) {
		ArrayList<ArrayList<ArrayList<Double>>> ans = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			ans.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				ans.get(state).add(new ArrayList<>());
				for(int residuePos=0;residuePos<objFcnSearch[state][subState].settings.mutRes.size();++residuePos) {
					ans.get(state).get(subState).add(0.0);
				}
				ans.get(state).get(subState).trimToSize();
			}
			ans.get(state).trimToSize();
		}
		ans.trimToSize();
		return ans;
	}

	private void clearPos2Score() {
		for(int state=0;state<pos2Score.size();++state) {
			for(int subState=0;subState<pos2Score.get(state).size();++subState) {
				for(int residuePos=0;residuePos<pos2Score.get(state).get(subState).size();++residuePos) {
					pos2Score.get(state).get(subState).set(residuePos, 0.0);
				}
			}
		}
	}

	private ArrayList<ArrayList<HashMap<Integer, HashSet<String>>>> allocatePos2AAs(MSSearchProblem[][] objFcnSearch) {
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

	private void clearPos2AAs() {
		for(int state=0;state<pos2AAs.size();++state) {
			for(int subState=0;subState<pos2AAs.get(state).size();++subState) {
				for(Integer residuePos : pos2AAs.get(state).get(subState).keySet()) {
					pos2AAs.get(state).get(subState).get(residuePos).clear();
				}
			}
		}
	}

	private boolean isMatchingResidue(Residue res, String resNum, String AAType) {
		return res.fullName.contains(resNum) && res.fullName.contains(AAType);
	}

	private boolean isMutable(Residue res, MSSearchProblem search) {
		return search.flexRes.contains(res.getPDBResNumber().trim());
	}

	protected boolean isUpperBoundPFProxy(BigDecimal objFcnCoeff, boolean isUnbound) {
		if(objFcnCoeff.compareTo(BigDecimal.ZERO)<0 && !isUnbound ||
				objFcnCoeff.compareTo(BigDecimal.ZERO)>0 && isUnbound)
			return true;
		return false;
	}

	protected void setPos2EnergyContrib(MSSearchProblem search, int state, int subState, double offset) {

		//double maxVal = search.settings.stericThreshold;
		double maxVal = 0;

		ArrayList<Integer> pos1s = search.getPosNums();
		for(int pos1 : pos1s) {
			double sumPairwise = 0;
			long numPairwise = 0;
			double pos1Val = 0;

			//also used pruned rcs at pos?
			ArrayList<Integer> pos1RCs = search.pruneMat.unprunedRCsAtPos(pos1);

			ArrayList<Integer> pos2s = pos1s;
			//ArrayList<Integer> pos2s = search.getPosNums();
			for(int pos2 : pos2s) {
				if(pos1==pos2) continue;

				ArrayList<Integer> pos2RCs = search.pruneMat.unprunedRCsAtPos(pos2);

				for(int rc1 : pos1RCs) {

					String rc1AAType = search.confSpace.posFlex.get(pos1).RCs.get(rc1).AAType;
					if(!pos2AAs.get(state).get(subState).get(pos1).contains(rc1AAType)) continue;

					for(int rc2 : pos2RCs) {

						String rc2AAType = search.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
						if(!pos2AAs.get(state).get(subState).get(pos2).contains(rc2AAType)) continue;

						//cap energy at steric threshold
						double pairwise = Math.min(search.emat.getPairwise(pos1, rc1, pos2, rc2), maxVal);
						sumPairwise += pairwise;
						numPairwise++;
					}
				}
			}

			pos1Val = sumPairwise/(double)numPairwise;
			if(Double.isNaN(pos1Val)) pos1Val = 0;

			pos2Score.get(state).get(subState).set(pos1, pos1Val);
		}
	}

	private void setPos2EnergyContrib(MSSearchProblem search, int state, 
			int subState, ScoredConf conf, MultiTermEnergyFunction mef, 
			double offset, boolean isUpperBoundPFProxy) {

		/*
		if(isUpperBoundPFProxy) {
			setPos2EnergyContrib(search, state, subState, offset);
			return;
		}
		*/

		ArrayList<Integer> unassignedPos = search.getPosNums(false);
		int[] assignments = conf.getAssignments();

		for(int pos : unassignedPos) {
			double intraAtPos = 0;
			double pairwiseAtPos = 0;

			// get intra contribution
			String resNum = search.flexRes.get(pos);//search.settings.mutRes.get(pos);
			int rcNum = assignments[pos];
			String AAType = search.confSpace.posFlex.get(pos).RCs.get(rcNum).AAType;

			for(int i=0;i<mef.getTerms().size();++i) {

				if( mef.getTerms().get(i) instanceof SingleResEnergy ) {
					Residue res = ((SingleResEnergy)mef.getTerms().get(i)).getRes();

					if(isMatchingResidue(res, resNum, AAType)) {
						intraAtPos += mef.getCoeffs().get(i) * mef.getTerms().get(i).getEnergy();
					}
				}

				else if( mef.getTerms().get(i) instanceof ResPairEnergy ) {
					Residue res1 = ((ResPairEnergy)mef.getTerms().get(i)).getRes1();
					Residue res2 = ((ResPairEnergy)mef.getTerms().get(i)).getRes2();

					//mutable to mutable residue interaction
					if( (isMatchingResidue(res1, resNum, AAType) && isMutable(res2, search)) ||
							(isMatchingResidue(res2, resNum, AAType) && isMutable(res1, search)) ) {
						double partialPairwise = mef.getCoeffs().get(i) * mef.getTerms().get(i).getEnergy();
						pairwiseAtPos += 0.5 * partialPairwise;
					}

					// See EnergyFunctionGenerator.intraAndShellEnergy(...)
					// intra energy consists of single residue energy and sum of pairs consisting of residue to shell energies
					else if( (isMatchingResidue(res1, resNum, AAType) && !isMutable(res2, search)) || 
							(isMatchingResidue(res2, resNum, AAType) && !isMutable(res1, search)) ) {						
						intraAtPos += mef.getCoeffs().get(i) * mef.getTerms().get(i).getEnergy();
					}
				}
			}

			double val = intraAtPos + pairwiseAtPos;
			pos2Score.get(state).get(subState).set(pos, val);
		}
	}

	ArrayList<ArrayList<String>> getAllowedAAs(int state, int subState) {
		ArrayList<ArrayList<String>> ans = new ArrayList<>();

		HashMap<Integer,HashSet<String>> p2AA = pos2AAs.get(state).get(subState);
		for(int pos=0;pos<p2AA.size();++pos) {
			ArrayList<String> tmp = new ArrayList<>();
			tmp.addAll(p2AA.get(pos));
			ans.add(tmp);
		}
		return ans;
	}

	private void computeStateScores(MSConfigFileParser cfp, BigDecimal objFcnCoeff, 
			MSSearchProblem search, int state, int subState, boolean isUnbound) {
		//find rigid gmec using conf enumeration

		//make conf search factory (i.e. A* tree)
		search.settings.overrideMultiSequence = true;
		ConfSearchFactory confSearchFactory = MSKStarFactory.makeConfSearchFactory(search, cfp);

		ArrayList<Integer> posNums = search.getPosNums(true);
		ArrayList<ArrayList<String>> AATypeOptions = getAllowedAAs(state, subState);

		ConfSearch tree = confSearchFactory.make(search.emat, search.getUpdatedPruningMatrix(posNums, AATypeOptions));
		search.settings.overrideMultiSequence = false;

		if(tree instanceof ConfAStarTree) ((ConfAStarTree)tree).stopProgress();
		ConfSearch.Splitter confsSplitter = new ConfSearch.Splitter(tree);
		ConfSearch.Splitter.Stream scoreConfs = confsSplitter.makeStream();

		//get rigid gmec
		ScoredConf conf = scoreConfs.next();
		if(conf == null) return;
		MultiTermEnergyFunction mef = search.getDecomposedEnergy(conf.getAssignments(), false);

		double offset = (mef.getPreCompE()-conf.getScore())/(double)search.getNumPos();
		/*
		if(Math.abs(offset)>0.1)
			throw new RuntimeException("ERROR: conf score "+ conf.getScore() +" != MEF energy " + mef.getPreCompE());
		 */

		boolean isUpperBoundPFProxy = isUpperBoundPFProxy(objFcnCoeff, isUnbound);
		setPos2EnergyContrib(search, state, subState, conf, mef, offset, isUpperBoundPFProxy);
	}

	private void computeStateScores(LMB objFcn, MSSearchProblem[][] objFcnSearch, KStarScore[] objFcnScores, boolean parallel) {		
		ArrayList<ResisueOrderWorker> workers = new ArrayList<>();
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				boolean isUnbound = subState != objFcnSearch[state].length-1;
				workers.add(new ResisueOrderWorker(state, subState, isUnbound));
			}
		}

		if(!parallel) {
			for(ResisueOrderWorker w : workers) computeStateScores(objFcnScores[w.state].getSettings().cfp, 
					objFcn.getCoeffs()[w.state], objFcnSearch[w.state][w.subState], 
					w.state, w.subState, w.isUnbound);
		}
		//execute in parallel
		else {
			workers.parallelStream().forEach(w -> computeStateScores(objFcnScores[w.state].getSettings().cfp, 
					objFcn.getCoeffs()[w.state], objFcnSearch[w.state][w.subState], 
					w.state, w.subState, w.isUnbound));
		}
	}

	protected void clearStoredStructures() {
		clearPos2Score();
		clearPos2AAs();
		ra2AAa.clear();
	}

	protected void storeResidue2AAAssignments(ResidueAssignment residueAssignment, 
			ArrayList<ArrayList<ArrayList<AAAssignment>>> aaAssignments) {
		ra2AAa.put(residueAssignment, aaAssignments);
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

	protected void computeAAsAndResAssignments(MSSearchProblem[][] objFcnSearch, ArrayList<Integer> unassignedPos, int numMaxMut) {
		//compute for assigned pos
		for(int state=0;state<objFcnSearch.length;++state) {
			for(int substate=0;substate<objFcnSearch[state].length;++substate) {
				for(int pos : objFcnSearch[state][substate].getPosNums(true)) {

					MSSearchProblem search = objFcnSearch[state][substate];

					if(search.settings.AATypeOptions.get(pos).size() != 1)
						throw new RuntimeException("ERROR: there must be exactly one AA per residue assignment");

					String allowedAA = search.settings.AATypeOptions.get(pos).get(0);

					pos2AAs.get(state).get(substate).get(pos).add(allowedAA);
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
							HashSet<String> aasAtPos = pos2AAs.get(state).get(substate).get(rPos);
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
						HashSet<String> aasAtPos = pos2AAs.get(state).get(substate).get(pos);
						if(aasAtPos.size() != 1)
							throw new RuntimeException("ERROR: there must be exactly one AA per residue assignment");
					}
				}
			}

			//2) each aa is contained in AA options in the search.settings
			for(int state=0;state<objFcnSearch.length;++state) {
				for(int substate=0;substate<objFcnSearch[state].length;++substate) {
					for(int pos : objFcnSearch[state][substate].getPosNums(false)) {
						HashSet<String> aasAtPos = pos2AAs.get(state).get(substate).get(pos);
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

	protected BigDecimal getStateScoreSum(int state, ResidueAssignment assignment) {		
		double ans = 0;
		//unbound states
		for(int subState=0;subState<assignment.length()-1;++subState) {
			ArrayList<Integer> unboundPos = assignment.get(subState);
			if(unboundPos.size()==0) continue;
			if(unboundPos.size()>1) throw new RuntimeException("ERROR: unbound state was split into more than one position");
			int pos = unboundPos.get(0);//contains at most one value		
			ans += pos2Score.get(state).get(subState).get(pos);
		}

		//bound state
		int subState = assignment.length()-1;
		ArrayList<Integer> boundPos = assignment.get(subState);
		for(int pos : boundPos) {
			ans += pos2Score.get(state).get(subState).get(pos);
		}

		return BigDecimal.valueOf(ans);
	}

	protected BigDecimal getStateScoreRatio(int state, ResidueAssignment assignment) {
		BigDecimal ans = BigDecimal.ONE.setScale(64, RoundingMode.HALF_UP);
		//unbound states
		for(int subState=0;subState<assignment.length()-1;++subState) {
			ArrayList<Integer> unboundPos = assignment.get(subState);
			if(unboundPos.size()==0) continue;
			if(unboundPos.size()>1) throw new RuntimeException("ERROR: unbound state was split into more than one position");
			int pos = unboundPos.get(0);//contains at most one value		
			ans = ans.multiply(boltzmann.calc(pos2Score.get(state).get(subState).get(pos)));
		}
		if(ans.compareTo(BigDecimal.ZERO)==0) ans = new BigDecimal("1e-63").setScale(64, RoundingMode.HALF_UP);

		//bound state
		int subState = assignment.length()-1;
		ArrayList<Integer> boundPos = assignment.get(subState);
		BigDecimal numerator = BigDecimal.ZERO.setScale(64, RoundingMode.HALF_UP);
		for(int pos : boundPos) {
			numerator = numerator.add(boltzmann.calc(pos2Score.get(state).get(subState).get(pos)));
		}

		ans = numerator.divide(ans, RoundingMode.HALF_UP);
		return ans;
	}


	protected BigDecimal getScore(ResidueAssignment assignment, 
			MSSearchProblem[][] objFcnSearch,
			BigDecimal[] coeffs,
			BigDecimal fScore) {

		int numStates = coeffs.length;
		BigDecimal ans = BigDecimal.ZERO.setScale(16, RoundingMode.HALF_UP);

		//get assignment score
		for(int state=0;state<numStates;++state) {
			//sign of 0 does not contribute to score
			if(coeffs[state].compareTo(BigDecimal.ZERO)==0) continue;
			//ans = ans.add(coeffs[state].multiply(getStateScoreRatio(state, assignment)));
			ans = ans.add(getStateScoreSum(state, assignment));
		}

		return ans;
	}

	protected ArrayList<ResidueAssignmentScore> scoreResidueAssignments(LMB objFcn, 
			MSSearchProblem[][] objFcnSearch,
			KStarScore[] objFcnScores) {

		BigDecimal[] stateVals = new BigDecimal[objFcnScores.length];
		for(int i=0;i<stateVals.length;++i) stateVals[i] = objFcnScores[i].getScore();
		BigDecimal fScore = objFcn.eval(stateVals);

		//now score each assignment
		ArrayList<ResidueAssignmentScore> assignmentScores = new ArrayList<>();
		for(ResidueAssignment residueAssignment : ra2AAa.keySet()) {
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
				return a1.score.compareTo(a2.score)<=0 ? -1 : 1;
			}
		});

		ResidueAssignment best = order.get(0).assignment;
		return best;
	}

	protected ArrayList<ArrayList<ArrayList<AAAssignment>>> getAAAssignments(ResidueAssignment residueAssignment) {
		ArrayList<ArrayList<ArrayList<AAAssignment>>> ans = ra2AAa.get(residueAssignment);
		return ans;
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

		computeStateScores(objFcn, objFcnSearch, objFcnScores, MSKStarNode.PARALLEL_EXPANSION);

		ArrayList<ResidueAssignmentScore> raScores = scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores);

		return raScores;
	}

	@Override
	public ArrayList<ArrayList<ArrayList<AAAssignment>>> getNextAssignments(LMB objFcn,
			MSSearchProblem[][] objFcnSearch, KStarScore[] objFcnScores, int numMaxMut) {

		ArrayList<ResidueAssignmentScore> raScores = scoreResidueAssignments(objFcn, objFcnSearch, objFcnScores, numMaxMut);

		ResidueAssignment best = getBestResidueAssignment(raScores);

		return getAAAssignments(best);
	}
}
