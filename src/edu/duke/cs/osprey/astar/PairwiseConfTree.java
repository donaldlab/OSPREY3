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
    protected void scoreNode(AStarNode node) {
    	
		EnergyMatrix emat = this.emat;
		int[][] unprunedRCsAtPos = this.unprunedRCsAtPos;
		
		int[] conf = node.getNodeAssignments();
		int[] undefinedRCIndices = node.getUndefinedRCIndices();
    
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
			int minRCIndex = -1;
			
			// for each conf...
			int[] rcs1 = unprunedRCsAtPos[pos1];
			int n1 = rcs1.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rcs1[j];
				
				double rcEnergy = getUndefinedRCEnergy(conf, pos1, rc1, j);
				
				if (rcEnergy < minRCEnergy) {
					minRCEnergy = rcEnergy;
					minRCIndex = j;
				}
			}
			
			hscore += minRCEnergy;
			undefinedRCIndices[pos1] = minRCIndex;
		}
		
    	resetSplitPositions();
    	
    	node.setScore(gscore + hscore);
    	
    	// DEBUG
    	assert (node.getScore() == scoreConf(node.getNodeAssignments()));
    }
    
    @Override
    protected void scoreNode(AStarNode parent, AStarNode node) {
    	// TODO
    	scoreNode(node);
    }
    
    //@Override
    protected double scoreConfDifferentialNope(AStarNode parent, int pos, int rc) {
    	
    	int[] parentConf = parent.getNodeAssignments();

    	double score = parent.getScore();
    	
    	// this should not be assigned in the parent
    	assert (parentConf[pos] < 0);
    	
    	// subtract the undefined contribution
    	int rcIndex = parent.getUndefinedRCIndices()[pos];
    	score -= getUndefinedRCEnergy(null, pos, unprunedRCsAtPos[pos][rcIndex], rcIndex);
    	
    	
    	// add the defined contribution
    	
    	// TEMP
		int[] conf = new int[numPos];
		System.arraycopy(parentConf, 0, conf, 0, numPos);
		conf[pos] = rc;
    	double expectedScore = scoreConf(conf);
    	double err = Math.abs(score - expectedScore);
    	err /= score;
    	assert (err <= 1e-15) : String.format("bad score: exp=%f, obs=%f", expectedScore, score);
    
    	return score;
    }
    
    @Override
    protected double getMinPairwiseEnergy(int[] conf, int pos1, int rc1, int rc1i, int pos2) {
		int index = this.precomputedMinOffsets[pos1*(pos1 - 1)/2 + pos2] + rc1i;
		return this.precomputedMins[index];
	}
}
