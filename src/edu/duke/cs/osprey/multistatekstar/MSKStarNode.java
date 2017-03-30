package edu.duke.cs.osprey.multistatekstar;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.multistatekstar.ResidueOrder.AAScore;
import edu.duke.cs.osprey.tools.ObjectIO;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Node for multistate k* tree
 *
 */
public class MSKStarNode {

	public static LMB OBJ_FUNC;//global objective function that we are optimizing for
	public static ArrayList<String[]> WT_SEQS;//wild-type bound state sequences for all states
	public static int NUM_MAX_MUT;//maximum number of allowed mutations

	private KStarScore[] ksLB;//lower bound k* objects
	private KStarScore[] ksUB;//upper bound k* objects
	private KStarScore[] ksObjFunc;//k* objects that minimize objective function
	private BigDecimal[] kss;//kstar score values
	private BigDecimal score;//objective function value; smaller is better

	public MSKStarNode(
			KStarScore[] ksLB, 
			KStarScore[] ksUB
			) {
		this.ksLB = ksLB;
		this.ksUB = ksUB;
		this.ksObjFunc = new KStarScore[ksLB.length];
		this.kss = new BigDecimal[ksLB.length];
		score = null;
	}

	public boolean isFullyAssigned() {
		return ksLB[0].isFullyAssigned();
	}

	public boolean isFinal() {
		return ksLB[0].isFinal();
	}

	private void computeKStarBounds(ArrayList<MSKStarNode> nodes, boolean parallel) {		
		if(!parallel) {
			KStarScore score;
			for(MSKStarNode node : nodes) {
				for(int state=0;state<ksLB.length;++state) {

					score = node.ksLB[state];
					if(score!=null) {
						score.compute(Integer.MAX_VALUE);
						if(!score.constrSatisfied()) break;
					}

					score = node.ksUB[state];
					if(score!=null) {
						score.compute(Integer.MAX_VALUE);
						if(!score.constrSatisfied()) break;
					}
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
		}
		nodes.removeAll(remove);
	}

	public ArrayList<MSKStarNode> split(ParamSet msParams) {
		ArrayList<MSKStarNode> ans = new ArrayList<>();
		if(isFullyAssigned()) {
			if(isFinal()) {

			}

			else {

			}
		}

		else {
			ksObjFunc = getStateKStarObjects(OBJ_FUNC);
			int numStates = ksObjFunc.length;
			MSSearchProblem[][] objFcnSearch = new MSSearchProblem[numStates][];
			for(int state=0;state<numStates;++state) objFcnSearch[state] = ksObjFunc[state].getSettings().search;	

			//residue ordering is determined by the objective function search problems, rather than just upper bound or lower bounds
			ResidueOrder order = ResidueOrderFactory.getResidueOrder(msParams, objFcnSearch);
			ArrayList<ArrayList<ArrayList<AAScore>>> splits = order.getNextAssignments(objFcnSearch, NUM_MAX_MUT-getNumMutations(0));

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
						break;//unbound state partition function upper bound(s) are 0
					}
					newKsLB[state] = lb;

					//make ub
					newKsUB[state] = split(ksUB[state], splits, splitIndex);
				}
				if(!addNode) continue;

				MSKStarNode child = new MSKStarNode(newKsLB, newKsUB);
				ans.add(child);
			}
		}

		computeKStarBounds(ans, false);
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
				pfs[subState] = parent.getPartitionFunction(subState);
			}
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

	private int getNumMutations(int state) {
		KStarScore score = ksLB[state];
		int boundState = score.getSettings().search.length-1;
		ArrayList<ArrayList<String>> AATypeOptions = score.getSettings().search[boundState].settings.AATypeOptions;

		int ans = 0;
		for(int pos : score.getSettings().search[boundState].getPosNums(true)) {
			if(AATypeOptions.get(pos).size()>0) ans++;
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
	
	public boolean constrSatisfiedGlobal() {
		for(KStarScore score : getStateKStarObjects(OBJ_FUNC)) {
			if(score==null) continue;
			if(!score.constrSatisfied()) return false;
		}
		return true;
	}

	public String toString() {
		return "";
	}
}
