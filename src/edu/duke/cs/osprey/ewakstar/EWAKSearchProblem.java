package edu.duke.cs.osprey.ewakstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.multistatekstar.UpdatePruningMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
public class EWAKSearchProblem extends SearchProblem implements UpdatePruningMatrix {

	private ArrayList<Integer> splitPosNums; 
	private ArrayList<ArrayList<String>> splitAAs;
	private PruningMatrix origPruneMat;

	public EWAKSearchProblem(SearchProblem other,
			ArrayList<Integer> splitPosNums, 
			ArrayList<ArrayList<String>> splitAAs) {
		super(other);
		this.origPruneMat = other.pruneMat;
		this.splitAAs = splitAAs;
		this.splitPosNums = splitPosNums;
	}

	public PruningMatrix updatePruningMatrix(
			ArrayList<Integer> splitPosNums, 
			ArrayList<ArrayList<String>> splitAAs
			) {
		UpdatedPruningMatrix ans = new UpdatedPruningMatrix(origPruneMat);
		for(int pos : splitPosNums) {
			for(int rc : origPruneMat.unprunedRCsAtPos(pos)) {
				String rcAAType = confSpace.posFlex.get(pos).RCs.get(rc).AAType;
				//is in split position, but not a desired AA type
				if(!splitAAs.get(pos).contains(rcAAType))
					ans.markAsPruned(new RCTuple(pos, rc));
			}
		}
		pruneMat = ans;

		return pruneMat;
	}

	public PruningMatrix updatePruningMatrix() {
		return updatePruningMatrix(splitPosNums, splitAAs);
	}

}
