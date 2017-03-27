package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

public class ResidueOrderStaticSequential implements ResidueOrder {

	protected ArrayList<ArrayList<PriorityQueue<AAScore>>> pqs;
	protected AAScore[][] nextPos;

	public ResidueOrderStaticSequential(MSSearchProblem[][] objFcnSearch) {
		pqs = null;
		nextPos = null;
		allocate(objFcnSearch);
		init(objFcnSearch);
	}
	
	private void allocate(MSSearchProblem[][] objFcnSearch) {
		nextPos = new AAScore[objFcnSearch.length][];
		for(int state=0;state<nextPos.length;++state) {
			nextPos[state] = new AAScore[objFcnSearch[state].length];
		}
		
		//create a priority queue for each state and substate
		for(int state=0;state<objFcnSearch.length;++state) {
			pqs.add(new ArrayList<>());
			for(int subState=0;subState<objFcnSearch[state].length;++subState)
				pqs.get(state).add(new PriorityQueue<AAScore>(14, new Comparator<AAScore>() {
					@Override
					public int compare(AAScore o1, AAScore o2) {
						return o1.score < o2.score ? -1 : 1;
					}		
				}));
			pqs.get(state).trimToSize();
		}
		pqs.trimToSize();
	}

	/**
	 * iterate through all aatypeoptions and assign an increasing numerical score
	 */
	protected void init(MSSearchProblem[][] objFcnSearch) {
		for(int state=0;state<pqs.size();++state) {
			int score = -1;
			for(int subState=0;subState<pqs.get(state).size();++subState) {
				ArrayList<ArrayList<String>> AATypeOptions = objFcnSearch[state][subState].settings.AATypeOptions;
				for(int residuePos=0;residuePos<AATypeOptions.size();++residuePos) {
					ArrayList<String> AATypes = AATypeOptions.get(residuePos);
					for(int AATypePos=0;AATypePos<AATypes.size();++AATypePos) {
						pqs.get(state).get(subState).add(new AAScore(residuePos, AATypePos, ++score));
					}
				}
			}
		}
	}

	@Override
	public AAScore[][] getNextAssignment(MSSearchProblem[][] objFcnSearch) {
		//substates with no assignable positions will have MAX_VALUE
		for(int state=0;state<nextPos.length;++state) {
			Arrays.fill(nextPos[state], null);
			for(int subState=0;subState<objFcnSearch[state].length;++subState) {
				MSSearchProblem search = objFcnSearch[state][subState];
				if(search.getNumUndefinedPos() > 0) {
					nextPos[state][subState] = pqs.get(state).get(subState).poll();
				}
			}
		}
		return nextPos;
	}
	
	protected int distanceToWT(String[] boundStateWT, MSSearchProblem boundState, AAScore aas) {
		int dist=0;
		for(int pos : boundState.getPosNums(true)) {
			if(!boundState.settings.AATypeOptions.get(pos).get(0).equalsIgnoreCase(boundStateWT[pos])) 
				dist++;
		}
		
		if(!boundState.allowedAAs.get(aas.residuePos).get(aas.AATypePos).equalsIgnoreCase(boundStateWT[aas.residuePos])) 
			dist++;
		
		return dist;
	}

}
