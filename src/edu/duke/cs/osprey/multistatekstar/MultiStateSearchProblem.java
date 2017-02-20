package edu.duke.cs.osprey.multistatekstar;

import edu.duke.cs.osprey.confspace.SearchProblem;

@SuppressWarnings("serial")
public class MultiStateSearchProblem extends SearchProblem {

	SearchProblemSettings spSet;
	
	public MultiStateSearchProblem(SearchProblem other) {
		super(other);
	}

	public boolean isFullyDefined() {
		return spSet==null || spSet.redFlexRes.size()==flexRes.size();
	}
	
	public MultiStateSearchProblem reduce(SearchProblemSettings spSet) {
		MultiStateSearchProblem other = new MultiStateSearchProblem(this);
		other.spSet = spSet;
		return other;
	}
}
