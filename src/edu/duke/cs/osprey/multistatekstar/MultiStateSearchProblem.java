package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;

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

	public void prunePmat(SearchProblem search, double pruningWindow, double stericThresh) {

		BigInteger numDesiredConfs = BigInteger.valueOf(65536);
		
		// single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = ((QPruningMatrix)pruneMat).countUpdates();
		int oldNumUpdates;

		Pruner dee = new Pruner(search, pruneMat, true, stericThresh, pruningWindow, 
				search.useEPIC, search.useTupExpForSearch);
		dee.setVerbose(false);

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = ((QPruningMatrix)pruneMat).countUpdates();
		} while (numUpdates > oldNumUpdates && 
				((QPruningMatrix)pruneMat).getNumReducedUnprunedConfs().compareTo(numDesiredConfs) > 0);
		
	}
}
