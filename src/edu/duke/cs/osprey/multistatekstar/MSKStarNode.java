package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfEnergyCalculator;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.multistatekstar.KStarScore.KStarScoreType;
import edu.duke.cs.osprey.multistatekstar.ResidueOrder.AAAssignment;
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

	public static ResidueOrder RESIDUE_ORDER;
	public static boolean PARALLEL_EXPANSION;
	public static boolean SUBLINEAR_AT_LEAF_NODES;
	public static boolean DEBUG = false;

	private KStarScore[] ksLB;//lower bound k* objects
	private KStarScore[] ksUB;//upper bound k* objects
	private KStarScore[] ksObjFunc;//k* objects that minimize objective function
	private BigDecimal[] kss;//kstar score values
	private BigDecimal score;//objective function value; smaller is better
	private int numPruned;
	public boolean isRoot;

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
		this.isRoot = false;
	}

	public String getSequence(int state) {
		int numSubStates = ksLB[state].getSettings().search.length;
		return ksLB[state].getSettings().search[numSubStates-1].settings.getFormattedSequence();
	}

	public int getNumStates() {
		return ksLB.length;
	}

	public int getNumSubStates() {
		return ksLB[0].getSettings().search.length;
	}

	public int getNumPruned() {
		return numPruned;
	}

	boolean keepUBScore() {
		KStarScore[] scores = getStateKStarObjects(OBJ_FUNC);
		int size = scores.length;

		for(int i=0; i<size; ++i) {
			KStarScore ofScore = scores[i];
			KStarScore lbScore = ksLB[i];
			if(ofScore != null && ofScore instanceof KStarScoreUpperBound && 
					ofScore.getDenom().compareTo(BigDecimal.ZERO)==0 &&
					lbScore.getDenom().compareTo(BigDecimal.ZERO)>0) {
				return true;
			}
		}
		return false;
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

	private void printDebugStatement(MSKStarNode node) {
		KStarScore[] parentScore = getStateKStarObjects(OBJ_FUNC);
		KStarScore[] childScore = node.getStateKStarObjects(OBJ_FUNC);
		
		for(int state=0;state<getNumStates();++state) {

			System.out.println();
			System.out.println("state: "+state);
			
			KStarScore parent = parentScore[state];
			System.out.println("parent score type: "+parent.getClass().getSimpleName());
			System.out.println("parent: "+parent.toString());
			
			KStarScore child = childScore[state];
			System.out.println("child score type: "+child.getClass().getSimpleName());
			System.out.println("child: "+child.toString());
			
			System.out.println();
		}
	}

	private boolean childScoreIsLessThanParentScore(MSKStarNode node) {
		//root node can have a null score, since all its children are valid
		if(this.getScore()==null) return false;
		
		if(node.getScore().compareTo(this.getScore())<0) {

			System.out.println(String.format("child-parent: %12e", node.getScore().subtract(this.getScore())));

			if(DEBUG) {
				printDebugStatement(node);
			}
			
			return true;
		}
		return false;
	}

	private void setChildScores(ArrayList<MSKStarNode> nodes, boolean parallel) {		
		if(!parallel) {
			KStarScore score;
			for(MSKStarNode node : nodes) {

				for(int state=0;state<ksLB.length;++state) {

					if(OBJ_FUNC.getCoeffs()[state].compareTo(BigDecimal.ZERO)==0)
						continue;

					score = node.ksLB[state];
					if(score!=null) {
						if(!score.isFinal()) {
							score.compute(Integer.MAX_VALUE);
						}
						else {
							score.computeUnboundStates(Integer.MAX_VALUE);
							score.computeBoundState(SUBLINEAR_AT_LEAF_NODES ? score.getSettings().ecalcs[0].getParallelism() : Integer.MAX_VALUE);
						}
					}

					score = node.ksUB[state];
					if(score!=null && !node.ksLB[state].equals(node.ksUB[state])) {
						if(!score.isFinal()) {
							score.compute(Integer.MAX_VALUE);
						}
						else {
							score.computeUnboundStates(Integer.MAX_VALUE);
							score.computeBoundState(SUBLINEAR_AT_LEAF_NODES ? score.getSettings().ecalcs[0].getParallelism() : Integer.MAX_VALUE);
						}
					}

				}
			}		
		} 

		else {
			ArrayList<KStarScore> scores = new ArrayList<>();
			for(MSKStarNode node : nodes) {
				for(int state=0;state<ksLB.length;++state) {

					if(OBJ_FUNC.getCoeffs()[state].compareTo(BigDecimal.ZERO)==0)
						continue;

					if(node.ksLB[state]!=null) {
						scores.add(node.ksLB[state]);
					}
					// for leaves, upper and lower bounds are the same, so don't 
					// try to compute both
					if(node.ksUB[state]!=null && !node.ksUB[state].equals(node.ksLB[state])) {
						scores.add(node.ksUB[state]);
					}
				}
			}
			scores.parallelStream().forEach(score -> {	
				if(!score.isFinal()) score.compute(Integer.MAX_VALUE);
				else{
					score.computeUnboundStates(Integer.MAX_VALUE);
					score.computeBoundState(SUBLINEAR_AT_LEAF_NODES ? score.getSettings().ecalcs[0].getParallelism() : Integer.MAX_VALUE);
				}
			});
		}

		//remove nodes that violate local constraints
		ArrayList<MSKStarNode> remove = new ArrayList<>();
		for(MSKStarNode node : nodes) {
			if(!node.constrSatisfiedLocal()) {
				remove.add(node);
				continue;
			}

			//this is an interesting corner case. if the upper bound score object
			//has a denom value of 0 (i.e. the partition function lower bounds are 0),
			//we must keep the node and expand it to a leaf. this is because the
			//lower bound partition function is defined on the rigid rotamer emat.
			//keeping the parent score ensures that we will process this node first
			else if(!this.isRoot && node.keepUBScore()) {
				node.setScore(this.getScore());
			}

			//set score based on k* bound scores
			else {
				node.setScore(OBJ_FUNC);
			}

			//check scores
			if(childScoreIsLessThanParentScore(node)) {
				throw new RuntimeException("ERROR: child score must be >= parent score");
			}
		}

		numPruned += remove.size();
		nodes.removeAll(remove);
	}

	public ArrayList<MSKStarNode> splitUnassigned() {
		if(isFullyAssigned()) throw new RuntimeException("ERROR: cannot split a "
				+ "fully assigned node unless continuous minimization is enabled");

		ArrayList<ArrayList<ArrayList<AAAssignment>>> splits = RESIDUE_ORDER.getNextAssignments(
				OBJ_FUNC, 
				getStateKStarSearch(OBJ_FUNC),
				getStateKStarObjects(OBJ_FUNC),
				NUM_MAX_MUT-getNumMutations(0)
				);

		ArrayList<MSKStarNode> ans = new ArrayList<>();	
		int numStates = getNumStates();
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
				//special case: if child is final, then ub = lb
				newKsUB[state] = lb.getSettings().isFinal ? lb : split(ksUB[state], splits, splitIndex);
			}
			if(!addNode) continue;

			MSKStarNode child = new MSKStarNode(newKsLB, newKsUB);
			ans.add(child);
		}

		setChildScores(ans, PARALLEL_EXPANSION);
		return ans;
	}

	private KStarScore split(KStarScore parent, ArrayList<ArrayList<ArrayList<AAAssignment>>> splits, int index) {
		//substate / split index / assignments for split
		MSKStarSettings kSet = new MSKStarSettings(parent.getSettings());

		int numSubStates = kSet.search.length;
		PartitionFunction[] pfs = new PartitionFunction[numSubStates];
		Arrays.fill(pfs, null);

		for(int subState=0;subState<numSubStates;++subState) {
			ArrayList<AAAssignment> assignments = splits.get(subState).get(index);
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

		//else, make a bounding type of k* object
		//in the minimized case, a fully assigned k* object is also a bounding
		//type of object
		else {
			return MSKStarFactory.makeKStarScore(kSet, pfs);
		}
	}

	private MSSearchProblem splitSearch(int subState, MSSearchProblem parent, ArrayList<AAAssignment> splits) {
		//make new search settings
		MSSearchSettings sSet = (MSSearchSettings) ObjectIO.deepCopy(parent.settings);

		for(AAAssignment aa : splits) {
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
					numPruned++;
					break;//unbound state partition function upper bound(s) are 0
				}
				//compute a tiny bit of the bound state
				//default to 1 * getparallelism
				lb.computeBoundState(SUBLINEAR_AT_LEAF_NODES ? lb.getSettings().ecalcs[0].getParallelism() : Integer.MAX_VALUE);
				newKsLB[state] = lb;

				//ub = lb
				newKsUB[state] = lb;
			}

			if(!addNode) return ans;
			child = new MSKStarNode(newKsLB, newKsUB);
		}

		else {
			// this applies to discrete AND continuous cases
			child = this;
			for(KStarScore lb : child.ksLB) {//when final, lb objects are the same as ub objects
				if(lb.isComputed()) continue;
				lb.computeBoundState(SUBLINEAR_AT_LEAF_NODES ? lb.getSettings().ecalcs[0].getParallelism() : Integer.MAX_VALUE);
			}
		}

		if(child.constrSatisfiedLocal()) {
			BigDecimal parentScore = this.getScore();

			child.setScore(OBJ_FUNC);

			if(DEBUG) {
				BigDecimal childScore = child.getScore();
				if(childScore.compareTo(parentScore)<0)
					throw new RuntimeException(String.format("ERROR: in leaf self-expansion, childScore: %12e < parentScore: %12e", childScore, parentScore));
			}

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

		ans.getSettings().computeGMEC = true;

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

	public void setScore(BigDecimal val) {
		score = val;
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

	public MSSearchProblem[][] getStateKStarSearch(LMB lmb) {
		KStarScore[] scores = getStateKStarObjects(lmb);
		int numStates = scores.length;
		MSSearchProblem[][] search = new MSSearchProblem[numStates][];
		for(int state=0;state<numStates;++state) search[state] = scores[state].getSettings().search;
		return search;
	}

	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Score: "+String.format("%12e", getScore())+"\n");
		KStarScore[] scores = getStateKStarObjects(OBJ_FUNC);
		for(int state=0;state<scores.length;++state) {
			KStarScore score = scores[state];
			sb.append("State"+state+": "+score.toString()+"\n");
		}
		return sb.toString();
	}
}
