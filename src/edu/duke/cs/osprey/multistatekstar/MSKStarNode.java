package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.ResidueOrder.AAScore;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	//static parameters are set by the tree object
	public static LMB OBJ_FUNC;//global objective function that we are optimizing for
	public static ArrayList<String[]> WT_SEQS;//wild-type bound state sequences for all states
	public static int NUM_MAX_MUT;//maximum number of allowed mutations
	public static ParamSet MS_PARAMS;//multistate parameters
	public static SearchProblem SEARCH_CONT[][];//continuous search problems
	public static SearchProblem SEARCH_DISC[][];//discrete search problems
	public static ConfEnergyCalculator.Async[][] ECALCS_CONT;//energy calculators for continuous emats
	public static ConfEnergyCalculator.Async[][] ECALCS_DISC;//energy calculators for discrete emats

	public static ResidueOrder RESIDUE_ORDER = null;
	private static int PARALLELISM_MULTIPLIER = 1;

	private KStarScore[] ksLB;//lower bound k* objects
	private KStarScore[] ksUB;//upper bound k* objects
	private KStarScore[] ksObjFunc;//k* objects that minimize objective function
	private BigDecimal[] kss;//kstar score values
	private BigDecimal score;//objective function value; smaller is better
	private int numPruned;

	public MSKStarNode(
			KStarScore[] ksLB, 
			KStarScore[] ksUB
			) {
		this.ksLB = ksLB;
		this.ksUB = ksUB;
		this.ksObjFunc = new KStarScore[ksLB.length];
		this.kss = new BigDecimal[ksLB.length];
		this.score = null;
		this.numPruned = 0;
	}

	public int getNumPruned() {
		return numPruned;
	}
	
	public int getNumAssignedResidues() {
		MSSearchProblem[] search = ksLB[0].getSettings().search;
		return search[search.length-1].getNumAssignedPos();
	}
	
	public boolean isFullyAssigned() {
		return ksLB[0].isFullyAssigned();
	}

	public boolean isFinal() {
		return ksLB[0].isFinal();
	}

	public boolean isLeafNode() {
		//kssLB is never null for fully defined sequences, minimized or not
		for(int i=0;i<ksLB.length;++i) if(!ksLB[i].isFullyProcessed()) return false;
		return true;
	}

	private boolean constrSatisfiedLocal() {
		for(int state=0;state<ksLB.length;++state) {
			if(!ksLB[state].constrSatisfied()) return false;
			if(!ksUB[state].constrSatisfied()) return false;
		}
		return true;
	}

	public boolean constrSatisfiedLocalObjFunc() {
		for(KStarScore score : getStateKStarObjects(OBJ_FUNC)) {
			if(score==null) continue;
			if(!score.constrSatisfied()) return false;
		}
		return true;
	}

	private void setNodeScores(ArrayList<MSKStarNode> nodes, boolean parallel) {
		if(!parallel) {
			KStarScore score;
			for(MSKStarNode node : nodes) {
				for(int state=0;state<ksLB.length;++state) {
					score = node.ksLB[state];
					if(score!=null) score.compute(Integer.MAX_VALUE);
					score = node.ksUB[state];
					if(score!=null) score.compute(Integer.MAX_VALUE);
				}
			}	
		} 

		else {
			ArrayList<KStarScore> scores = new ArrayList<>();
			for(MSKStarNode node : nodes) {
				for(int state=0;state<ksLB.length;++state) {
					if(node.ksLB[state]!=null) scores.add(node.ksLB[state]);
					if(node.ksUB[state]!=null) scores.add(node.ksUB[state]);
				}
			}
			scores.parallelStream().forEach(score -> score.compute(Integer.MAX_VALUE));
		}

		//remove nodes that violate local constraints
		ArrayList<MSKStarNode> remove = new ArrayList<>();
		for(MSKStarNode node : nodes) {
			if(!node.constrSatisfiedLocal()) remove.add(node);
			else node.setScore(OBJ_FUNC);
		}
		numPruned += remove.size();
		nodes.removeAll(remove);
	}

	public ArrayList<MSKStarNode> splitUnassigned() {
		if(isFullyAssigned()) throw new RuntimeException("ERROR: cannot split a "
				+ "fully assigned node unless continuous minimization is enabled");

		ArrayList<MSKStarNode> ans = new ArrayList<>();

		ksObjFunc = getStateKStarObjects(OBJ_FUNC);
		int numStates = ksObjFunc.length;
		MSSearchProblem[][] objFcnSearch = new MSSearchProblem[numStates][];
		for(int state=0;state<numStates;++state) objFcnSearch[state] = ksObjFunc[state].getSettings().search;	

		//residue ordering is determined by the objective function search problems, rather than just upper bound or lower bounds
		RESIDUE_ORDER = RESIDUE_ORDER == null ? ResidueOrderFactory.getResidueOrder(MS_PARAMS, objFcnSearch) : RESIDUE_ORDER;
		ArrayList<ArrayList<ArrayList<AAScore>>> splits = RESIDUE_ORDER.getNextAssignments(objFcnSearch, NUM_MAX_MUT-getNumMutations(0));

		//each split is applied to every state ub and lb
		int numSplits = splits.get(0).size();
		for(int splitIndex=0;splitIndex<numSplits;++splitIndex) {
			KStarScore[] newKsLB = new KStarScore[numStates];
			KStarScore[] newKsUB = new KStarScore[numStates];

			boolean addNode = true;
			//for each split, make a new child node
			for(int state=0;state<numStates;++state) {
				//make lb
				KStarScore lb = split(ksLB[state], splits, splitIndex);
				lb.computeUnboundStates(Integer.MAX_VALUE);
				if(!lb.constrSatisfied()) {
					addNode = false;
					numPruned++;
					break;//unbound state partition function upper bound(s) are 0
				}
				newKsLB[state] = lb;

				//make ub
				//special case: if discrete and child is final, ub = lb
				newKsUB[state] = lb.getSettings().isFinal ? lb : split(ksUB[state], splits, splitIndex);
			}
			if(!addNode) continue;

			MSKStarNode child = new MSKStarNode(newKsLB, newKsUB);
			ans.add(child);
		}

		setNodeScores(ans, true);
		return ans;
	}

	private KStarScore split(KStarScore parent, ArrayList<ArrayList<ArrayList<AAScore>>> splits, int index) {
		//substate / split index / assignments for split
		MSKStarSettings kSet = new MSKStarSettings(parent.getSettings());
		
		int numSubStates=kSet.search.length;
		PartitionFunction[] pfs = new PartitionFunction[numSubStates];
		Arrays.fill(pfs, null);

		for(int subState=0;subState<numSubStates;++subState) {
			ArrayList<AAScore> assignments = splits.get(subState).get(index);
			if(assignments.size()>0)//update search problem,otherwise we keep parent search problem
				kSet.search[subState] = splitSearch(subState, kSet.search[subState], assignments);
			else {
				if(subState==numSubStates-1) throw new RuntimeException("ERROR: there must always be a split for the bound state");
				//no split for substate, so re-use partition function
				pfs[subState] = parent.getPartitionFunction(subState);
			}
		}

		//special case: if discrete and parent is not fully assigned but child is, 
		//then create a non-bounding type of k* object
		ParamSet sParams = kSet.cfp.getParams();
		if(!sParams.getBool("DOMINIMIZE") && kSet.search[numSubStates-1].isFullyAssigned()) {
			//new search problems, keeping newly created search problem settings
			MSSearchProblem[] search = new MSSearchProblem[numSubStates];
			for(int subState=0;subState<numSubStates;++subState) {
				MSSearchSettings sSet = kSet.search[subState].settings;
				search[subState] = new MSSearchProblem(SEARCH_DISC[kSet.state][subState], sSet);
			}
			
			return MSKStarFactory.makeKStarScore(
					MS_PARAMS, kSet.state, kSet.cfp, kSet.constraints,
					null, search,
					null, ECALCS_DISC[kSet.state], KStarScoreType.Discrete
					);
		}

		return MSKStarFactory.makeKStarScore(kSet, pfs);
	}

	private MSSearchProblem splitSearch(int subState, MSSearchProblem parent, ArrayList<AAScore> splits) {
		//make new search settings
		MSSearchSettings sSet = (MSSearchSettings) ObjectIO.deepCopy(parent.settings);
		for(AAScore aa : splits) {
			//update mutres
			sSet.mutRes.set(aa.residuePos, parent.flexRes.get(aa.residuePos));
			//restrict allowed AAs
			String AAType = sSet.AATypeOptions.get(aa.residuePos).get(aa.AATypePos);
			sSet.AATypeOptions.get(aa.residuePos).clear();
			sSet.AATypeOptions.get(aa.residuePos).add(AAType);
			sSet.AATypeOptions.get(aa.residuePos).trimToSize();
		}

		MSSearchProblem ans = new MSSearchProblem(parent, sSet);
		return ans;
	}

	public ArrayList<MSKStarNode> splitFullyAssigned() {
		//applies only to minimized confs
		if(!ksLB[0].getSettings().cfp.getParams().getBool("DOMINIMIZE"))
			throw new RuntimeException("ERROR: can only split a fully assigned node when continuous minimization is enabled");

		MSKStarNode child = null;
		ArrayList<MSKStarNode> ans = new ArrayList<>();
		if(!isFinal()) {
			//transition from k* lower and upper bound to k* with continuous min
			int numStates = ksLB.length;
			KStarScore[] newKsLB = new KStarScore[numStates];
			KStarScore[] newKsUB = new KStarScore[numStates];

			boolean addNode = true;
			for(int state=0;state<numStates;++state) {
				//make lb
				KStarScore lb = splitToMinimized(ksLB[state]);
				lb.computeUnboundStates(Integer.MAX_VALUE);
				if(!lb.constrSatisfied()) {
					addNode = false;
					break;//unbound state partition function upper bound(s) are 0
				}
				//compute a tiny bit of the bound state
				//default to 16 * getparallelism
				lb.computeBoundState(PARALLELISM_MULTIPLIER * lb.getSettings().ecalcs[0].getParallelism());
				newKsLB[state] = lb;

				//ub = lb
				newKsUB[state] = lb;
			}

			if(!addNode) return ans;
			child = new MSKStarNode(newKsLB, newKsUB);
		}

		else {
			child = this;
			for(KStarScore lb : child.ksLB) {
				if(lb.isComputed()) continue;
				lb.computeBoundState(PARALLELISM_MULTIPLIER * lb.getSettings().ecalcs[0].getParallelism());
			}
		}

		if(child.constrSatisfiedLocal()) {
			child.setScore(OBJ_FUNC);
			ans.add(child);
		}

		return ans;
	}

	/**
	 * fully assigned parent creates a fully assigned child
	 * @param parent
	 * @return
	 */
	private KStarScore splitToMinimized(KStarScore parent) {
		int state = parent.getSettings().state;
		int numSubStates = parent.getSettings().search.length;
		MSSearchProblem[] search = new MSSearchProblem[numSubStates];
		for(int subState=0;subState<numSubStates;++subState) {
			//same sequence as before
			MSSearchSettings sSet = parent.getSettings().search[subState].settings;
			search[subState] = new MSSearchProblem(SEARCH_CONT[state][subState], sSet);
		}

		KStarScore ans = MSKStarFactory.makeKStarScore(
				MS_PARAMS, state, parent.getSettings().cfp, parent.getSettings().constraints,
				search, null,
				ECALCS_CONT[state], null, KStarScoreType.Minimized
				);

		return ans;
	}

	private int getNumMutations(int state) {
		KStarScore score = ksLB[state];
		int boundState = score.getSettings().search.length-1;
		ArrayList<ArrayList<String>> AATypeOptions = score.getSettings().search[boundState].settings.AATypeOptions;

		int ans = 0;
		for(int pos : score.getSettings().search[boundState].getPosNums(true)) {
			if(AATypeOptions.get(pos).size() != 1) throw new RuntimeException("ERROR: assigned positions must have 1 AA");
			else if(!AATypeOptions.get(pos).get(0).equalsIgnoreCase(WT_SEQS.get(state)[pos]))
				ans++;
		}
		return ans;
	}

	public BigDecimal getScore() {
		return score;
	}

	public void setScore(LMB lmb) {
		score = lmb.eval(getStateKStarScores(lmb));
	}

	public BigDecimal[] getStateKStarScores(LMB lmb) {
		BigDecimal[] coeffs = lmb.getCoeffs();
		//always want a lower bound on the lmb value
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) kss[i] = ksUB[i].getUpperBoundScore();
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) kss[i] = ksLB[i].getLowerBoundScore();
			else {//coeffs[i]==0
				if(lmb.equals(OBJ_FUNC)) throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
				else kss[i] = BigDecimal.ZERO;
			}
		}
		return kss;
	}

	public KStarScore[] getStateKStarObjects(LMB lmb) {
		BigDecimal[] coeffs = lmb.getCoeffs();
		//always want a lower bound on the lmb value
		for(int i=0;i<coeffs.length;++i) {
			if(coeffs[i].compareTo(BigDecimal.ZERO) < 0) ksObjFunc[i] = ksUB[i];
			else if(coeffs[i].compareTo(BigDecimal.ZERO) > 0) ksObjFunc[i] = ksLB[i];
			else {//coeffs[i]==0
				if(lmb.equals(OBJ_FUNC)) throw new RuntimeException("ERROR: objective function coefficient cannot be 0");
				else ksObjFunc = null;
			}
		}
		return ksObjFunc;
	}

	public String toString() {
		KStarScore[] scores = getStateKStarObjects(OBJ_FUNC);
		StringBuilder sb = new StringBuilder();
		for(int state=0;state<scores.length;++state) {
			KStarScore score = scores[state];
			sb.append("State"+state+": "+score.toString()+"\n");
		}
		return sb.toString();
	}
}
