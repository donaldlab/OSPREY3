package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;

@SuppressWarnings("serial")
public class MultiStateSearchProblem extends SearchProblem {

	SearchSettings settings;

	public MultiStateSearchProblem(SearchProblem other, 
			SearchSettings spSet) {
		super(other);
		this.settings = spSet;
		this.pruneMat = getReducedPruningMatrix();
		this.allowedAAs = spSet.AATypeOptions;
		this.flexRes = spSet.mutRes;
	}
	
	public boolean isFullyDefined() {
		return settings==null || settings.mutRes.size()==confSpace.numPos;
	}

	private QPruningMatrix getReducedPruningMatrix() {
		return new QPruningMatrix(this, settings.mutRes, settings.AATypeOptions);
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
