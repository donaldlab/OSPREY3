package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

public class ResidueOrder {

	public class ResidueScore {
		public int residuePos;
		public int AATypePos;
		public double score;

		public ResidueScore(int residuePos, int AATypePos, double score) {
			this.residuePos = residuePos;
			this.AATypePos = AATypePos;
			this.score = score;
		}
	}

	public enum ResidueOrderType {
		StaticSequential,
		StaticMinDom,
		StaticObjFuncHMean,
		DynamicObjFuncHMean;
	}

	public class ResidueScoreComparator implements Comparator<ResidueScore> {
		@Override
		public int compare(ResidueScore o1, ResidueScore o2) {
			return o1.score < o2.score ? -1 : 1;
		}
	}

	private MSSearchProblem searchCont[][];
	private MSSearchProblem searchDisc[][];
	private ArrayList<ArrayList<PriorityQueue<ResidueScore>>> scoreQueues;

	public ResidueOrder(
			MSSearchProblem searchCont[][],
			MSSearchProblem searchDisc[][]
			) {
		this.searchCont = searchCont;
		this.searchDisc = searchDisc;
		MSSearchProblem[][] tmp = this.searchCont != null ? this.searchCont : this.searchDisc;

		//create a priority queue for each state and substate
		Comparator<ResidueScore> comparator = new ResidueScoreComparator();
		for(int state=0;state<tmp.length;++state) {
			scoreQueues.add(new ArrayList<>());
			for(int subState=0;subState<tmp[state].length;++subState)
				scoreQueues.get(state).add(new PriorityQueue<ResidueScore>(comparator));
			scoreQueues.get(state).trimToSize();
		}
		scoreQueues.trimToSize();
	}

	public void scoreResidues(ResidueOrderType type) {
		switch(type) {
		case StaticSequential:
			scoreStaticSequential();
			break;
		case StaticMinDom:
		case StaticObjFuncHMean:
		case DynamicObjFuncHMean:
		default:
			throw new UnsupportedOperationException("ERROR: unsupported type "+type);
		}
	}
	
	/**
	 * iterate through all aatypeoptions and assign an increasing numerical score
	 */
	private void scoreStaticSequential() {
		MSSearchProblem[][] tmp = this.searchCont != null ? this.searchCont : this.searchDisc;
		int score=-1;
		for(int state=0;state<scoreQueues.size();++state) {
			for(int subState=0;subState<scoreQueues.get(state).size();++subState) {
				ArrayList<ArrayList<String>> AATypeOptions = tmp[state][subState].settings.AATypeOptions;
				for(int residuePos=0;residuePos<AATypeOptions.size();++residuePos) {
					ArrayList<String> AATypes = AATypeOptions.get(residuePos);
					for(int AATypePos=0;AATypePos<AATypes.size();++AATypePos) {
						scoreQueues.get(state).get(subState).add(new ResidueScore(residuePos, AATypePos, ++score));
					}
				}
			}
		}
	}
	
	/**
	 * get the next residue and AA assignment for the specified state and substate
	 * @param state
	 * @param subState
	 * @return
	 */
	public ResidueScore getNextAssignment(int state, int subState) {
		ResidueScore ans = scoreQueues.get(state).get(subState).poll();
		//make sure position is not already assigned
		if(searchCont!=null && searchCont[state][subState].flexRes.get(ans.residuePos).equals("-1"))
			throw new RuntimeException("ERROR: residue position "+ans.residuePos+" is already assigned");
		if(searchDisc!=null && searchDisc[state][subState].flexRes.get(ans.residuePos).equals("-1"))
			throw new RuntimeException("ERROR: residue position "+ans.residuePos+" is already assigned");
		return ans;
	}

}
