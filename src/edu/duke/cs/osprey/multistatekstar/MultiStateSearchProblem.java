package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.confspace.SearchProblem;

@SuppressWarnings("serial")
public class MultiStateSearchProblem extends SearchProblem {

	SearchProblemSettings spSet;
	
	public MultiStateSearchProblem(SearchProblem other, 
			SearchProblemSettings spSet) {
		super(other);
		this.spSet = spSet;
		this.pruneMat = getReducedPruningMatrix();
	}

	public boolean isFullyDefined() {
		return spSet==null || spSet.mutRes.size()==flexRes.size();
	}
	
	private QPruningMatrix getReducedPruningMatrix() {
		return new QPruningMatrix(this, spSet.mutRes, spSet.AATypeOptions);
	}
}
