package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class TraditionalPairwiseHScorer implements AStarScorer {
	
	private EnergyMatrix emat;
	private RCs rcs;
	
	private double[][][] undefinedEnergies; // indexed by (pos1,pos2), rc at pos1
	private ConfAStarNode cachedNode;
	private double[][] cachedEnergies;
	
	public TraditionalPairwiseHScorer(EnergyMatrix emat, RCs rcs) {
		this.emat = emat;
		this.rcs = rcs;
		
		int numPos = emat.getNumPos();
		
		// pre-compute all undefined energy terms
		undefinedEnergies = new double[numPos][][];
		for (int pos1=0; pos1<numPos; pos1++) {
			
			int numRCs = rcs.get(pos1).length;
			undefinedEnergies[pos1] = new double[numRCs][];
			
			for (int i=0; i<numRCs; i++) {
				int rc1 = rcs.get(pos1)[i];
				
				undefinedEnergies[pos1][i] = new double[numPos];
			
				for (int pos2=0; pos2<pos1; pos2++) {
					
					// compute the min over rc2
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						minEnergy = Math.min(minEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					undefinedEnergies[pos1][i][pos2] = minEnergy;
				}
			}
		}
		
		// allocate space for the cache
		cachedNode = null;
		cachedEnergies = new double[numPos][];
		for (int pos=0; pos<numPos; pos++) {
			cachedEnergies[pos] = new double[rcs.get(pos).length];
		}
	}
	
	public TraditionalPairwiseHScorer make() {
		return new TraditionalPairwiseHScorer(emat, rcs);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
    	// bound energy of undefined conf
    	double hscore = 0;
    	
    	calcCachedEnergies(confIndex, rcs);
    	
		// for each undefined pos...
		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos = confIndex.undefinedPos[i];
			
			// find the lowest-energy rc at this pos
			double minRCEnergy = Double.POSITIVE_INFINITY;
			for (int j=0; j<rcs.get(pos).length; j++) {
				minRCEnergy = Math.min(minRCEnergy, cachedEnergies[pos][j]);
			}
			
			hscore += minRCEnergy;
		}
		
		return hscore;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations can make an impact
		
		// if the intermediate energies aren't cached, calculate them now
		if (cachedNode != confIndex.node) {
			calcCachedEnergies(confIndex, rcs);
			cachedNode = confIndex.node;
		}
		
    	// compute the h-score
    	double hscore = 0;
    	for (int i=0; i<confIndex.numUndefined; i++) {
    		int pos = confIndex.undefinedPos[i];
    		
    		// don't score at nextPos, it's defined now
    		if (pos == nextPos) {
    			continue;
    		}
    		
    		// compute the new min energy over all rcs
    		double minRCEnergy = Double.POSITIVE_INFINITY;
    		
    		double[] cachedEnergiesAtPos = cachedEnergies[pos];
    		double[][] undefinedEnergiesAtPos = undefinedEnergies[pos];
    		
			// for each rc at this pos...
			int[] rcsAtPos = rcs.get(pos);
			int n = rcsAtPos.length;
			for (int j=0; j<n; j++) {
				int rc = rcsAtPos[j];
				
				double rcEnergy = cachedEnergiesAtPos[j];
				
				// subtract undefined contribution
				if (pos > nextPos) {
					rcEnergy -= undefinedEnergiesAtPos[j][nextPos];
				}
				
				// add defined contribution
				rcEnergy += emat.getPairwise(pos, rc, nextPos, nextRc);
				
				minRCEnergy = Math.min(minRCEnergy, rcEnergy);
			}
			
			hscore += minRCEnergy;
    	}
    	
    	return hscore;
	}

	private void calcCachedEnergies(ConfIndex confIndex, RCs rcs) {
		
		// for each undefined pos...
		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos1 = confIndex.undefinedPos[i];
			
			// for each rc...
			int[] rcs1 = rcs.get(pos1);
			int n1 = rcs1.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rcs1[j];
				
				// start with the one-body energy
				double energy = emat.getOneBody(pos1, rc1);
				
				// add defined energies
				for (int k=0; k<confIndex.numDefined; k++) {
					int pos2 = confIndex.definedPos[k];
					int rc2 = confIndex.definedRCs[k];
					
					energy += emat.getPairwise(pos1, rc1, pos2, rc2);
				}
				
				// add undefined energies
				double[] energies = undefinedEnergies[pos1][j];
				for (int k=0; k<confIndex.numUndefined; k++) {
					int pos2 = confIndex.undefinedPos[k];
					if (pos2 < pos1) {
						energy += energies[pos2];
					}
				}
				
				cachedEnergies[pos1][j] = energy;
			}
		}
	}
}
