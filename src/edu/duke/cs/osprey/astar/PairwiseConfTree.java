package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

// this conf tree sacrifices higher-order tuples for much much faster speed =)
public class PairwiseConfTree extends ConfTree<FullAStarNode> { // TODO: SlimAStarNode

	private double[][][] undefinedEnergies; // indexed by (pos1,pos2), rc at pos1
	
	private AStarNode cachedNode;
	private double[][] cachedEnergies;
    
	public PairwiseConfTree(SearchProblem sp){
		this(sp, sp.pruneMat, sp.useEPIC);
	}

	public PairwiseConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC) {
		super(new FullAStarNode.Factory(sp.confSpace.numPos), sp, pruneMat, useEPIC);
		
		if (emat.hasHigherOrderTerms()) {
			throw new Error("Don't use PairwiseConfTree with higher order energy terms");
		}
		
		// pre-compute all undefined energy terms
		undefinedEnergies = new double[numPos][][];
		for (int pos1=0; pos1<numPos; pos1++) {
			
			int numRCs = unprunedRCsAtPos[pos1].length;
			undefinedEnergies[pos1] = new double[numRCs][];
			
			for (int i=0; i<numRCs; i++) {
				int rc1 = unprunedRCsAtPos[pos1][i];
				
				undefinedEnergies[pos1][i] = new double[numPos];
			
				for (int pos2=0; pos2<pos1; pos2++) {
					
					// compute the min over rc2
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rc2 : unprunedRCsAtPos[pos2]) {
						minEnergy = Math.min(minEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					undefinedEnergies[pos1][i][pos2] = minEnergy;
				}
			}
		}
		
		// allocate cache space
		cachedNode = null;
		cachedEnergies = new double[numPos][];
		for (int pos=0; pos<numPos; pos++) {
			cachedEnergies[pos] = new double[unprunedRCsAtPos[pos].length];
		}
	}
	
    @Override
    protected double scoreNode(FullAStarNode node) {
    	
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
    	
    	computeCachedEnergies(node);
    	
		// for each undefined pos...
		for (int i=0; i<numUndefined; i++) {
			int pos = undefinedPos[i];
			
			// find the lowest-energy rc at this pos
			double minRCEnergy = Double.POSITIVE_INFINITY;
			for (int j=0; j<unprunedRCsAtPos[pos].length; j++) {
				minRCEnergy = Math.min(minRCEnergy, cachedEnergies[pos][j]);
			}
			
			hscore += minRCEnergy;
		}
		
    	resetSplitPositions();
    	
    	node.setGScore(gscore);
    	node.setScore(gscore + hscore);
    	
    	// DEBUG
    	//assertScore(scoreConf(node.getNodeAssignments()), node.getScore());
    	
    	return node.getScore();
    }
    
    @Override
    protected double scoreNodeDifferential(FullAStarNode parent, FullAStarNode child, int nextPos, int nextRc) {
    	
    	// NOTE: child can be null, eg for scoreExpansionLevel()
    	
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
    	// however, copying this.xxx references to the stack empirically didn't help this time
    	
    	// update cached info
    	if (parent != cachedNode) { // yes, compare by reference
    		
    		int[] parentConf = parent.getNodeAssignments();
    		splitPositions(parentConf);
    		computeCachedEnergies(parent);
    		cachedNode = parent;
    		
    		// this should not be assigned in the parent
    		assert (parentConf[nextPos] < 0);
    	}
    	
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
    		
    		// compute the new min energy over all rcs
    		double minRCEnergy = Double.POSITIVE_INFINITY;
    		
    		double[] cachedEnergiesAtPos = cachedEnergies[pos];
    		double[][] undefinedEnergiesAtPos = undefinedEnergies[pos];
    		
			// for each rc at this pos...
			int[] rcs = unprunedRCsAtPos[pos];
			int n = rcs.length;
			for (int j=0; j<n; j++) {
				int rc = rcs[j];
				
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
    	
    	double score = gscore + hscore;
    	
    	if (child != null) {
			child.setGScore(gscore);
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
    protected double scoreConfDifferential(FullAStarNode parent, int nextPos, int nextRc) {
    	return scoreNodeDifferential(parent, null, nextPos, nextRc);
    }
    
    private void computeCachedEnergies(FullAStarNode node) {
    	assertSplitPositions();
    	
		// for each undefined pos...
		for (int i=0; i<numUndefined; i++) {
			int pos1 = undefinedPos[i];
			
			// for each rc...
			int[] rcs1 = unprunedRCsAtPos[pos1];
			int n1 = rcs1.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rcs1[j];
				
				// start with the one-body energy
				double energy = emat.getOneBody(pos1, rc1);
				
				// add defined energies
				for (int k=0; k<numDefined; k++) {
					int pos2 = definedPos[k];
					int rc2 = definedRCs[k];
					
					energy += emat.getPairwise(pos1, rc1, pos2, rc2);
				}
				
				// add undefined energies
				double[] energies = undefinedEnergies[pos1][j];
				for (int k=0; k<numUndefined; k++) {
					int pos2 = undefinedPos[k];
					if (pos2 < pos1) {
						energy += energies[pos2];
					}
				}
				
				cachedEnergies[pos1][j] = energy;
			}
		}
    }
    
    @SuppressWarnings("unused")
	private void assertScore(double expected, double observed) {
    	double err = Math.abs(expected - observed);
    	if (observed != 0) {
    		err /= observed;
    	}
    	assert (err <= 1e-12) : String.format("bad score:\nexp=%.18f\nobs=%.18f\nerr=%.18f", expected, observed, err);
    }
}
