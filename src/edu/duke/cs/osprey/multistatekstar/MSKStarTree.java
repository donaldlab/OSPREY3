package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.ResidueOrder.ResidueOrderType;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * multi state k* tree
 *
 */
public class MSKStarTree {

	int numTreeLevels;//number of residues with sequence
	//changes+1 level if we are doing continuous minimization

	LMB objFcn;//we are minimizing objFcn
	LMB[] msConstr;
	LMB[][] sConstr;

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

	PriorityQueue<MSKStarNode> pq;

	int numSeqsWanted;
	int numSeqsReturned;
	
	int numExpanded;
	int numPruned;

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

		this.numTreeLevels = numTreeLevels;
		this.objFcn = objFcn;
		this.msConstr = msConstr;
		this.sConstr = sConstr;
		this.AATypeOptions = AATypeOptions;
		this.numMaxMut = numMaxMut;
		this.numSeqsWanted = numSeqsWanted;
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
		
		numExpanded = 0;
		numPruned = 0;
		numSeqsReturned = 0;
		pq = null;
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
		if(curNode.isLeafNode()) {
			//state-specific constraints already checked
		} else {
			//check state-specific constraints
		}
		
		//check global constraints
		return false;
	}
	
	private ArrayList<MSKStarNode> getChildren(MSKStarNode curNode) {
		//for each state and substate, pick next position to expand

		//create search problems

		//compute and score children
		ArrayList<MSKStarNode> ans = new ArrayList<>();
		ans.trimToSize();
		return ans;
	}

	private MSKStarNode rootNode() {
		//decide whether to make upper or lower bound parition functions based on
		//the required objective function and constraints
		boolean[] makeLB = new boolean[numStates]; Arrays.fill(makeLB, false);
		boolean[] makeUB = new boolean[numStates]; Arrays.fill(makeUB, false);
		
		for(int state=0;state<numStates;++state) {
			//decide whether to create ub or lb according to objFcn
			//need lb, so make lb
			if(objFcn.getCoeffs()[state].compareTo(BigDecimal.ZERO) < 0) makeLB[state] = true;
			//need ub so make ub
			else if(objFcn.getCoeffs()[state].compareTo(BigDecimal.ZERO) > 0) makeUB[state] = true;
			
			//decide whether to make ub or lb based on global constraints
			for(LMB constr : msConstr) {
				//want to eliminate by lb, so make ub
				if(constr.getCoeffs()[state].compareTo(BigDecimal.ZERO) < 0) makeUB[state] = true;
				//want to eliminate by ub, so make lb
				else if(constr.getCoeffs()[state].compareTo(BigDecimal.ZERO) > 0) makeLB[state] = true;
			}
			
			//decide whether to make ub or lb based on state-specific constraints
			for(LMB constr : sConstr[state]) {
				int numSubStates = sConstr[state].length;
				for(int subState=0; subState<numSubStates; ++subState) {
					
					if(subState <= numSubStates-2) {
						//want to eliminate unbound states by lb, so make lb score, which has upper bound unbound state pfs
						if(constr.getCoeffs()[subState].compareTo(BigDecimal.ZERO) < 0) makeLB[state] = true;
						//want to eliminate unbound states by ub, so make ub score, which has lower bound unbound state pfs
						if(constr.getCoeffs()[subState].compareTo(BigDecimal.ZERO) > 0) makeUB[state] = true;
					}
					
					else if(subState == numSubStates-1) {//bound state partition function
						//want to eliminate bound state by lb, so make ub score, which has upper bound bound state pf
						if(constr.getCoeffs()[subState].compareTo(BigDecimal.ZERO) < 0) makeUB[state] = true;
						//want to eliminate bound state by ub, so make lb score, which has lower bound bound state pf
						else if(constr.getCoeffs()[subState].compareTo(BigDecimal.ZERO) > 0) makeLB[state] = true;
					}
				}
			}
		}

		KStarScore[] kssLB = new KStarScore[numStates]; Arrays.fill(kssLB, null);
		KStarScore[] kssUB = new KStarScore[numStates]; Arrays.fill(kssUB, null);
		KStarScoreType[] types = null;
		
		for(int state=0;state<numStates;++state) {
			boolean doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");
			if(doMinimize)
				types = new KStarScoreType[]{KStarScoreType.MinimizedLowerBound, KStarScoreType.MinimizedUpperBound};
			else
				types = new KStarScoreType[]{KStarScoreType.DiscreteLowerBound, KStarScoreType.DiscreteUpperBound};

			//adjust types by the kinds of partition functions that we want to create
			if(!makeLB[state]) types[0] = null;
			if(!makeUB[state]) types[1] = null;
			
			KStarScore[] scores = getRootKStarScores(state, types);
			kssLB[state] = scores[0];
			kssUB[state] = scores[1];
		}

		MSKStarNode ans = new MSKStarNode(this, kssLB, kssUB);
		ans.setScore(MSKStarNode.NEGATIVE_INFINITY);
		return ans;
	}

	private KStarScore[] getRootKStarScores(int state, KStarScoreType[] types) {
		boolean doMinimize = cfps[state].getParams().getBool("DOMINIMIZE");
		//[0] is lb, [1] is ub
		KStarScore[] ans = new KStarScore[types.length]; Arrays.fill(ans, null);

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

	public String[] nextSeq() {
		
		if(pq==null)
			initQueue(rootNode());
	
		MSKStarNode curNode;
		while(true) {
			curNode = pq.poll();
			
			if(curNode==null) {
				System.out.println("Multi-State K* tree empty...returning empty signal");
				return null;
			}
		}
	}

}
