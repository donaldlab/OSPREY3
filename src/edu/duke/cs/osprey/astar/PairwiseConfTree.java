package edu.duke.cs.osprey.astar;

import java.util.Arrays;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

// this conf tree sacrifices higher-order tuples for much much faster speed =)
public class PairwiseConfTree extends ConfTree {

    int[] precomputedMinOffsets;
    double[] precomputedMins;
    
	public PairwiseConfTree(SearchProblem sp){
		this(sp, sp.pruneMat, sp.useEPIC);
	}

	public PairwiseConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC) {
		super(sp, pruneMat, useEPIC);
		
		if (emat.hasHigherOrderTerms()) {
			throw new Error("Don't use PairwiseConfTree with higher order energy terms");
		}
    	
        // precompute some mins to speed up scoring
        precomputedMinOffsets = null;
        precomputedMins = null;
		int numPos = emat.getNumPos();
		int index = 0;
		int offset = 0;
		precomputedMinOffsets = new int[numPos*(numPos - 1)/2];
		Arrays.fill(precomputedMinOffsets, 0);
		for (int res1=0; res1<numPos; res1++) {
			for (int res2=0; res2<res1; res2++) {
				precomputedMinOffsets[index++] = offset;
				offset += unprunedRCsAtPos[res1].length;
			}
		}
		precomputedMins = new double[offset];
		Arrays.fill(precomputedMins, Double.POSITIVE_INFINITY);
		index = 0;
		for (int res1=0; res1<numPos; res1++) {
			for (int res2=0; res2<res1; res2++) {
				for (int rc1 : unprunedRCsAtPos[res1]) {
					
					// compute the min
					double energy = Double.POSITIVE_INFINITY;
					for (int rc2 : unprunedRCsAtPos[res2]) {
						energy = Math.min(energy, emat.getPairwise(res1, rc1, res2, rc2));
					}
					
					precomputedMins[index] = energy;
					index++;
				}
			}
		}
    }
    
    @Override
    protected double scoreNode(AStarNode node) {
    	
    	// OPTIMIZATION: this function is really only called once for the root node
    	// no need to optimize it
    	
		EnergyMatrix emat = this.emat;
		int[][] unprunedRCsAtPos = this.unprunedRCsAtPos;
		
		int[] conf = node.getNodeAssignments();
    
    	splitPositions(conf);
    	int numDefined = this.numDefined;
    	int numUndefined = this.numUndefined;
    	int[] definedPos = this.definedPos;
    	int[] definedRCs = this.definedRCs;
		int[] undefinedPos = this.undefinedPos;
		RCTuple rcTuple = this.rcTuple;
		
		// compute energy of defined conf
    	double gscore = emat.getConstTerm();
    	
    	// one body energies
		for (int i=0; i<numDefined; i++) {
			int pos1 = definedPos[i];
			int rc1 = definedRCs[i];
			
			gscore += emat.getOneBody(pos1, rc1);
		}
		
		// pairwise energies
		assert (numDefined == rcTuple.size());
		for (int i=0; i<numDefined; i++) {
			int pos1 = definedPos[i];
			int rc1 = definedRCs[i];
		
			for (int j=0; j<i; j++) {
				int pos2 = definedPos[j];
				int rc2 = definedRCs[j];
				
				gscore += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
    	// then bound energy of undefined conf
    	double hscore = 0;
    	
		// for each undefined pos...
		for (int i=0; i<numUndefined; i++) {
			int pos1 = undefinedPos[i];
			
			// find the lowest-energy rc at this pos
			double minRCEnergy = Double.POSITIVE_INFINITY;
			
			int[] rcs1 = unprunedRCsAtPos[pos1];
			int n1 = rcs1.length;
			double[] rcEnergies = new double[n1];
			node.setUndefinedRCEnergies(pos1, rcEnergies);
			
			// for each rc...
			for (int j=0; j<n1; j++) {
				int rc1 = rcs1[j];
				
				double rcEnergy = getUndefinedRCEnergy(conf, pos1, rc1, j);
				rcEnergies[j] = rcEnergy;
				
				minRCEnergy = Math.min(minRCEnergy, rcEnergy);
			}
			
			hscore += minRCEnergy;
		}
		
    	resetSplitPositions();
    	
    	node.setGScore(gscore);
    	node.setHScore(hscore);
    	node.setScore(gscore + hscore);
    	
    	// DEBUG
    	//assertScore(scoreConf(node.getNodeAssignments()), node.getScore());
    	
    	return node.getScore();
    }
    
    @Override
    protected double scoreNodeDifferential(AStarNode parent, AStarNode child, int nextPos, int nextRc) {
    	
    	// NOTE: child can be null, eg for scoreExpansionLevel()
    	
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
    	// however, copying this.xxx references to the stack empirically didn't help this time
		
    	int[] parentConf = parent.getNodeAssignments();
    	splitPositions(parentConf);
    	
    	// this should not be assigned in the parent
    	assert (parentConf[nextPos] < 0);
    	
    	// update the g-score
    	double gscore = parent.getGScore();
    	
    	// add the new one-body energy
    	gscore += emat.getOneBody(nextPos, nextRc);
    	
    	// add the new pairwise energies
    	for (int i=0; i<numDefined; i++) {
    		int pos = definedPos[i];
    		int rc = definedRCs[i];
    		gscore += emat.getPairwise(pos, rc, nextPos, nextRc);
    	}
    	
    	// compute the h-score
    	double hscore = 0;
    	for (int i=0; i<numUndefined; i++) {
    		int pos = undefinedPos[i];
    		
    		// don't score at nextPos, it's defined now
    		if (pos == nextPos) {
    			continue;
    		}
    		
    		// get partially-computed energies from the parent node
    		double[] undefinedRCEnergies = parent.getUndefinedRCEnergies(pos);
    		
    		// and copy them to the child
    		if (child != null) {
    			undefinedRCEnergies = undefinedRCEnergies.clone();
    			child.setUndefinedRCEnergies(pos, undefinedRCEnergies);
    		}
    		
    		// compute the new min energy over all rcs
    		double minRCEnergy = Double.POSITIVE_INFINITY;
    		
			// for each rc at this pos...
			int[] rcs = unprunedRCsAtPos[pos];
			int n = rcs.length;
			for (int j=0; j<n; j++) {
				int rc = rcs[j];
				
				double rcEnergy = undefinedRCEnergies[j];
				
				// subtract undefined contribution
				if (pos > nextPos) {
					rcEnergy -= getMinPairwiseEnergy(null, pos, rc, j, nextPos);
				}
				
				// add defined contribution
				rcEnergy += emat.getPairwise(pos, rc, nextPos, nextRc);
				
				// save updated energies to child
				if (child != null) {
					undefinedRCEnergies[j] = rcEnergy;
				}
				
				minRCEnergy = Math.min(minRCEnergy, rcEnergy);
			}
			
			hscore += minRCEnergy;
    	}
    	
    	resetSplitPositions();
    	
    	double score = gscore + hscore;
    	
    	if (child != null) {
			child.setGScore(gscore);
			child.setHScore(hscore);
			child.setScore(score);
    	}
    	
    	/* DEBUG
		int[] conf = new int[numPos];
		System.arraycopy(parentConf, 0, conf, 0, numPos);
		conf[newPos] = newRc;
		assertScore(scoreConf(conf), score);
		*/
    
    	return score;
    }
    
    @Override
    protected double scoreConfDifferential(AStarNode parent, int nextPos, int nextRc) {
    	return scoreNodeDifferential(parent, null, nextPos, nextRc);
    }
    
    @Override
    protected double getMinPairwiseEnergy(int[] conf, int pos1, int rc1, int rc1i, int pos2) {
    	assert (pos2 < pos1);
		int index = this.precomputedMinOffsets[pos1*(pos1 - 1)/2 + pos2] + rc1i;
		return this.precomputedMins[index];
	}
    
    @SuppressWarnings("unused")
	private void assertScore(double expected, double observed) {
    	double err = Math.abs(expected - observed);
    	err /= observed;
    	assert (err <= 1e-13) : String.format("bad score:\nexp=%.18f\nobs=%.18f\nerr=%.18f", expected, observed, err);
    }
}
