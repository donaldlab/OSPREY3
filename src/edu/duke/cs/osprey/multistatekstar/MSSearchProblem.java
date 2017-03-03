package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMethod;

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

	private QPruningMatrix prunePmat(SearchProblem search, double pruningWindow, double stericThresh) {

		QPruningMatrix qpm = (QPruningMatrix) pruneMat;
		//don't want to overprune
		BigInteger numDesiredConfs = BigInteger.valueOf(65536);
		
		//single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = qpm.countUpdates();
		int oldNumUpdates;

		Pruner dee = new Pruner(search, qpm, true, stericThresh, pruningWindow, 
				search.useEPIC, search.useTupExpForSearch);
		dee.setVerbose(false);

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			//pairs pruning can be time consuming for larger systems
			if(dee.enumerateCandidates(PruningMethod.getMethod("GOLDSTEIN PAIRS FULL")).size() < 32768)
				dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = qpm.countUpdates();
			
		} while (numUpdates > oldNumUpdates && 
				qpm.getNumReducedUnprunedConfs().compareTo(numDesiredConfs) > 0);
		
		return qpm;
	}
}
