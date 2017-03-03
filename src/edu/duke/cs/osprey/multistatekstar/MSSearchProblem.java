package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;

@SuppressWarnings("serial")
public class MSSearchProblem extends SearchProblem {

	SearchSettings settings;

	public MSSearchProblem(SearchProblem other, 
			SearchSettings settings) {
		super(other);
		if(settings==null) throw new RuntimeException("ERROR: search settings cannot be null");
		this.settings = settings;
		this.pruneMat = getReducedPruningMatrix();
		//this.allowedAAs = settings.AATypeOptions;
		//this.flexRes = settings.mutRes;
		prunePmat(this, settings.stericThreshold, settings.stericThreshold);
	}
	
	public boolean isFullyDefined() {
		return settings.mutRes.size()==confSpace.numPos;
	}

	private QPruningMatrix getReducedPruningMatrix() {
		return new QPruningMatrix(this, settings.mutRes, settings.AATypeOptions);
	}

	private void prunePmat(SearchProblem search, double pruningWindow, double stericThresh) {

		BigInteger numDesiredConfs = BigInteger.valueOf(65536);
		
		//single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = ((QPruningMatrix)pruneMat).countUpdates();
		int oldNumUpdates;

		Pruner dee = new Pruner(search, pruneMat, true, stericThresh, pruningWindow, 
				search.useEPIC, search.useTupExpForSearch);
		dee.setVerbose(false);
		
		//limit amount of time (in minutes) spent pruning
		double start = System.currentTimeMillis(), duration, timeout = 60;

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = ((QPruningMatrix)pruneMat).countUpdates();
			
			duration = (System.currentTimeMillis()-start)/(60*1000.0);
			
		} while (numUpdates > oldNumUpdates && 
				((QPruningMatrix)pruneMat).getNumReducedUnprunedConfs().compareTo(numDesiredConfs) > 0 &&
				duration < timeout);
		
	}
}
