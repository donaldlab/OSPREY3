package edu.duke.cs.osprey.multistatekstar;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Stopwatch;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * multi state k* tree
 *
 */
public class MSKStarTree {

	LMB objFcn;//we are minimizing objFcn
	LMB[] msConstr;
	LMB[][] sConstr;

	ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums;
	//mutable2StateResNum.get(state) maps levels in this tree to flexible positions for state

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;
	// MultiStateKStarTreeNode.assignments Assigns each level an index in 
	// AATypeOptions.get(level), and thus an AA type
	//If -1, then no assignment yet

	int numTreeLevels;//maximum number of mutable/flexible residues.
	//+1 if minimization is allowed
	
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

	PriorityQueue<MSKStarNode> pq;

	int numSeqsWanted;
	int numSeqsReturned;

	int numExpanded;
	int numSelfExpanded;
	int numPruned;
	
	Stopwatch stopwatch;

	public MSKStarTree(
			int numTreeLevels,
			int numStates,
			int numMaxMut,
			int numSeqsWanted,
			LMB objFcn,
			LMB[] msConstr,
			LMB[][] sConstr,
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

		this.objFcn = objFcn;
		this.msConstr = msConstr;
		this.sConstr = sConstr;
		this.AATypeOptions = AATypeOptions;
		this.numMaxMut = numMaxMut;
		this.numSeqsWanted = numSeqsWanted;
		this.wtSeqs = wtSeqs;
		this.numTreeLevels = numTreeLevels;
		this.numStates = numStates;
		this.searchCont = searchCont;
		this.searchDisc = searchDisc;
		this.ecalcsCont = ecalcsCont;
		this.ecalcsDisc = ecalcsDisc;
		this.mutable2StateResNums = mutable2StateResNums;

		this.cfps = cfps;
		this.msParams = msParams;

		this.numExpanded = 0;
		this.numSelfExpanded = 0;
		this.numPruned = 0;
		this.numSeqsReturned = 0;
		this.pq = null;
		
		this.stopwatch = new Stopwatch().start();
	}

	private void initQueue(MSKStarNode node) {
		pq = new PriorityQueue<MSKStarNode>(1024, new Comparator<MSKStarNode>() {
			@Override
			public int compare(MSKStarNode m1, MSKStarNode m2) {
				return m1.getScore().compareTo(m2.getScore()) < 0 ? -1 : 1;
			}
		});

		pq.add(node);
	}

	private boolean canPrune(MSKStarNode curNode) {
		//first check whether local constraints are satisfied
		//if(!curNode.constrSatisfiedLocalObjFunc()) return true;

		//now check global constraints
		for(LMB lmb : msConstr) {
			if(lmb.eval(curNode.getStateKStarScores(lmb)).compareTo(BigDecimal.ZERO) > 0)
				return true;
		}
		return false;
	}

	private ArrayList<MSKStarNode> getChildren(MSKStarNode curNode) {
		ArrayList<MSKStarNode> ans = new ArrayList<>();

		//pick next position to expand
		if(!curNode.isFullyAssigned())
			ans.addAll(curNode.splitUnassigned());

		else
			ans.addAll(curNode.splitFullyAssigned());

		ans.trimToSize();
		return ans;
	}

	private void initNodeStaticVars() {
		//initialize MSKStarNode
		MSKStarNode.OBJ_FUNC = this.objFcn;
		MSKStarNode.WT_SEQS = this.wtSeqs;
		MSKStarNode.NUM_MAX_MUT = this.numMaxMut;
		MSKStarNode.MS_PARAMS = this.msParams;
		MSKStarNode.SEARCH_CONT = this.searchCont;
		MSKStarNode.SEARCH_DISC = this.searchDisc;
		MSKStarNode.ECALCS_CONT = this.ecalcsCont;
		MSKStarNode.ECALCS_DISC = this.ecalcsDisc;
		MSKStarNode.RESIDUE_ORDER = null;
	}

	private MSKStarNode getRootNode() {

		initNodeStaticVars();

		KStarScore[] kssLB = new KStarScore[numStates];
		KStarScore[] kssUB = new KStarScore[numStates];
		KStarScoreType[] types = null;

		for(int state=0;state<numStates;++state) {
			boolean doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");
			if(doMinimize)
				types = new KStarScoreType[]{KStarScoreType.MinimizedLowerBound, KStarScoreType.MinimizedUpperBound};
			else
				types = new KStarScoreType[]{KStarScoreType.DiscreteLowerBound, KStarScoreType.DiscreteUpperBound};

			KStarScore[] scores = getRootKStarBounds(state, types);
			kssLB[state] = scores[0];
			kssUB[state] = scores[1];
		}

		MSKStarNode ans = new MSKStarNode(kssLB, kssUB);
		ans.setScore(objFcn);//set score per the objective function
		return ans;
	}

	private KStarScore[] getRootKStarBounds(int state, KStarScoreType[] types) {
		boolean doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");
		//[0] is lb, [1] is ub
		KStarScore[] ans = new KStarScore[types.length];

		ParamSet sParams = cfps[state].getParams();
		int numPartFuncs = sParams.getInt("NUMUBSTATES")+1;

		for(int i=0;i<types.length;++i) {

			KStarScoreType type = types[i];
			if(type == null) continue;

			MSSearchProblem[] seqSearchCont = doMinimize ? new MSSearchProblem[numPartFuncs] : null;
			MSSearchProblem[] seqSearchDisc = new MSSearchProblem[numPartFuncs];

			for(int subState=0;subState<numPartFuncs;++subState) {

				MSSearchSettings spSet = new MSSearchSettings();
				spSet.AATypeOptions = AATypeOptions.get(state).get(subState);
				ArrayList<String> mutRes = new ArrayList<>();
				for(int j=0;j<mutable2StateResNums.get(state).get(subState).size();++j) mutRes.add("-1");
				mutRes.trimToSize();
				spSet.mutRes = mutRes;
				spSet.stericThreshold = sParams.getDouble("STERICTHRESH");
				spSet.pruningWindow = sParams.getDouble("IVAL") + sParams.getDouble("EW");

				if(doMinimize) seqSearchCont[subState] = new MSSearchProblem(searchCont[state][subState], spSet);
				seqSearchDisc[subState] = new MSSearchProblem(searchDisc[state][subState], (MSSearchSettings) ObjectIO.deepCopy(spSet));
			}

			ans[i] = MSKStarFactory.makeKStarScore(
					msParams, state, cfps[state], sConstr[state],
					seqSearchCont, seqSearchDisc,
					ecalcsCont[state], ecalcsDisc[state], type
					);
		}

		return ans;
	}

	public String nextSeq() {

		if(pq==null)
			initQueue(getRootNode());

		MSKStarNode curNode;
		while(true) {
			curNode = pq.poll();

			if(curNode==null) {
				System.out.println("Multi-State K* tree empty...returning empty signal");
				return null;
			}

			//count number pruned by local constraints
			numPruned += curNode.getNumPruned();
			
			if(canPrune(curNode)) {
				numPruned++;
				continue;
			}

			else {		
				if(curNode.isLeafNode()) {
					numSeqsReturned++;
					reportProgress(curNode);
					return curNode.toString();
				}
				
				//expand
				ArrayList<MSKStarNode> children = getChildren(curNode);
				//expansion is either a refinement of the same node or creation
				//of completely new nodes
				if(children.size()==0 && curNode.equals(children.get(0))) numSelfExpanded++;
				numExpanded++;

				pq.addAll(children);
				
				//if(numExpanded % 8==0) 
					reportProgress(curNode);
			}
		}
	}
	
	private void reportProgress(MSKStarNode curNode) {
		MemoryUsage heapMem = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		System.out.println(String.format("level: %d/%d, score: %12e, expanded: %5d, pruned: %5d, seqs: %d/%d, time: %6s, heapMem: %.0f%%",
				curNode.getNumAssignedResidues(), numTreeLevels, curNode.getScore(), numExpanded, numPruned,
				numSeqsReturned, numSeqsWanted, stopwatch.getTime(2), 100f*heapMem.getUsed()/heapMem.getMax()));
	}

}
