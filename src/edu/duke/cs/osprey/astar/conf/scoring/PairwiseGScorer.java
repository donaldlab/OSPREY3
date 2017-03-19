package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class PairwiseGScorer implements AStarScorer {
	
	private EnergyMatrix emat;
	
	public PairwiseGScorer(EnergyMatrix emat) {
		this.emat = emat;
	}
	
	@Override
	public PairwiseGScorer make() {
		return new PairwiseGScorer(emat);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// constant term
    	double gscore = emat.getConstTerm();
    	
    	// one body energies
		for (int i=0; i<confIndex.getNumDefined(); i++) {
			int pos1 = confIndex.getDefinedPos()[i];
			int rc1 = confIndex.getDefinedRCs()[i];
			
			gscore += emat.getOneBody(pos1, rc1);
		}
		
		// pairwise energies
		for (int i=0; i<confIndex.getNumDefined(); i++) {
			int pos1 = confIndex.getDefinedPos()[i];
			int rc1 = confIndex.getDefinedRCs()[i];
		
			for (int j=0; j<i; j++) {
				int pos2 = confIndex.getDefinedPos()[j];
				int rc2 = confIndex.getDefinedRCs()[j];
				
				gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
		return gscore;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
    	// modify the parent node's g-score
    	double gscore = confIndex.getNode().getGScore();
    	
    	// add the new one-body energy
    	gscore += emat.getOneBody(nextPos, nextRc);
    	
    	// add the new pairwise energies
    	for (int i=0; i<confIndex.getNumDefined(); i++) {
    		int pos = confIndex.getDefinedPos()[i];
    		int rc = confIndex.getDefinedRCs()[i];
    		gscore += emat.getPairwise(pos, rc, nextPos, nextRc);
    	}
    	
    	return gscore;
	}
}
