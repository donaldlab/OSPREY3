package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class NaiveTraditionalPairwiseHScorer implements AStarScorer {
	
	private EnergyMatrix emat;
	
	public NaiveTraditionalPairwiseHScorer(EnergyMatrix emat) {
		this.emat = emat;
	}
	
	@Override
	public NaiveTraditionalPairwiseHScorer make() {
		return new NaiveTraditionalPairwiseHScorer(emat);
	}
	
	@Override
	public double calc(ConfIndex index, RCs rcs) {
		
		double hscore = 0;
		
		// get the score for each undefined position
		for (int i=0; i<index.numUndefined; i++) {
			int pos1 = index.undefinedPos[i];
			
			// min over possible assignments to pos1
			double pos1Score = Double.POSITIVE_INFINITY;
			for (int rc1 : rcs.get(pos1)) {
				
				double rcContrib = emat.getOneBody(pos1, rc1);

				// interactions with defined residues
				for (int j=0; j<index.numDefined; j++) {
					int pos2 = index.definedPos[j];
					int rc2 = index.definedRCs[j];
					rcContrib += emat.getPairwise(pos1, rc1, pos2, rc2);
				}

				// interactions with undefined residues
				for (int j=0; j<index.numUndefined; j++) {
					int pos2 = index.undefinedPos[j];
					if (pos2 >= pos1) {
						break;
					}

					// min over possible assignments to pos2
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
						minEnergy = Math.min(minEnergy, pairwiseEnergy);
					}

					rcContrib += minEnergy;
				}

				pos1Score = Math.min(pos1Score, rcContrib);
			}
		
			hscore += pos1Score;
		}
		
		return hscore;
	}
}
