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
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];
			
			gscore += emat.getOneBody(pos1, rc1);
		}
		
		// pairwise energies
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos1 = confIndex.definedPos[i];
			int rc1 = confIndex.definedRCs[i];
		
			for (int j=0; j<i; j++) {
				int pos2 = confIndex.definedPos[j];
				int rc2 = confIndex.definedRCs[j];
				
				gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
		return gscore;
	}

	public double calc(int[] assignments) {

		// constant term
		double gscore = emat.getConstTerm();

		for (int pos1=0; pos1<assignments.length; pos1++) {
			int rc1 = assignments[pos1];

			// one body energy
			gscore += emat.getOneBody(pos1, rc1);

			// pairwise energies
			for (int pos2=0; pos2<pos1; pos2++) {
				int rc2 = assignments[pos2];

				gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}

		return gscore;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
    	// modify the parent node's g-score
    	double gscore = confIndex.node.getGScore();
    	
    	// add the new one-body energy
    	gscore += emat.getOneBody(nextPos, nextRc);
    	
    	// add the new pairwise energies
    	for (int i=0; i<confIndex.numDefined; i++) {
    		int pos = confIndex.definedPos[i];
    		int rc = confIndex.definedRCs[i];
    		gscore += emat.getPairwise(pos, rc, nextPos, nextRc);
    	}
    	
    	return gscore;
	}
}
