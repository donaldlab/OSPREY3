package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class UpperBoundScorer implements AStarScorer{
	
	private EnergyMatrix emat;
	
	public UpperBoundScorer(EnergyMatrix emat) {
		this.emat = emat;
	}
	
	@Override
	public UpperBoundScorer make() {
		return new UpperBoundScorer(emat);
	}
	
	@Override
	public double calc(ConfIndex index, RCs rcs) {
		double uscore = 0;
		
		// get the score for each undefined position
				for (int i=0; i<index.numUndefined; i++) {
					int pos1 = index.undefinedPos[i];
					
					// max over possible assignments to pos1
					double pos1Score = Double.NEGATIVE_INFINITY;
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

							// max over possible assignments to pos2
							double maxEnergy = Double.NEGATIVE_INFINITY;
							for (int rc2 : rcs.get(pos2)) {
								double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
								maxEnergy = Math.max(maxEnergy, pairwiseEnergy);
							}

							rcContrib += maxEnergy;
						}

						pos1Score = Math.max(pos1Score, rcContrib);
					}
				
					uscore += pos1Score;
				}
		return uscore;
	}
}
