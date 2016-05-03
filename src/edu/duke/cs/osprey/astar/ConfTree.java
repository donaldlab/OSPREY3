/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author mhall44
 */
public class ConfTree extends AStarTree {
    //This implementation of an A* tree is intended for conformational search
    //AStarNode.nextAssignment is an array of length numPos; each position
    //stores the assigned RC, or -1 to indicate an unassigned position
    //this class supports both "traditional" (static, simple heuristic) A*
    //and improvements like dynamic A*
    //we may also want to allow other negative indices, to indicate partially assigned RCs

    int numPos;
    EnergyMatrix emat;
    
    
    int maxNumRCs;
    int[] unprunedRCsAtPos;
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue
    
    
    
    
    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    boolean traditionalScore = true;
    boolean useRefinement = false;//refine nodes (might want EPIC, MPLP, or something else)
    
    boolean useDynamicAStar = true;

    
    EPICMatrix epicMat = null;//to use in refinement
    ConfSpace confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined
    
    
    public ConfTree(SearchProblem sp){
        this(sp, sp.pruneMat, sp.useEPIC);
    }
    
    public ConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC){
        numPos = sp.confSpace.numPos;
        
        //see which RCs are unpruned and thus available for consideration
        
        // pack them into an efficient int matrix
        // we'll allocate a little too much space,
        // but we're talking numResidues*max(numRCs) here,
        // so it shouldn't be too big compared to the beastly tuple matrices
        
        // get the max number of RCs
        maxNumRCs = 0;
        for (int pos=0; pos<numPos; pos++) {
        	maxNumRCs = Math.max(maxNumRCs, pruneMat.unprunedRCsAtPos(pos).size());
        }
        assert (maxNumRCs > 0);
        
        // copy the RC numbers over
        unprunedRCsAtPos = new int[maxNumRCs*numPos];
        Arrays.fill(unprunedRCsAtPos, -1);
        for (int pos=0; pos<numPos; pos++) {
        	int index = pos*maxNumRCs;
        	for (int rc : pruneMat.unprunedRCsAtPos(pos)) {
        		assert (rc >= 0);
        		unprunedRCsAtPos[index++] = rc;
        	}
        }
        
        //get the appropriate energy matrix to use in this A* search
        if(sp.useTupExpForSearch)
            emat = sp.tupExpEMat;
        else {
            emat = sp.emat;
            
            if(useEPIC){//include EPIC in the search
                useRefinement = true;
                epicMat = sp.epicMat;
                confSpace = sp.confSpace;
                minPartialConfs = sp.epicSettings.minPartialConfs;
            }
        }
    }
    
    public int getNumRCs(int pos) {
    	int index = pos*maxNumRCs;
    	int count = 0;
    	for (int i=0; i<maxNumRCs; i++) {
    		if (unprunedRCsAtPos[index + i] < 0) {
    			break;
    		}
    		count++;
    	}
    	return count;
    }
    
    public BigInteger getNumConformations() {
    	BigInteger num = BigInteger.valueOf(1);
    	for (int pos=0; pos<numPos; pos++) {
    		num = num.multiply(BigInteger.valueOf(getNumRCs(pos)));
    	}
    	return num;
    }
    
    
    
    @Override
    public ArrayList<AStarNode> getChildren(AStarNode curNode) {
        
        if(isFullyAssigned(curNode))
            throw new RuntimeException("ERROR: Can't expand a fully assigned A* node");
        
        if(curNode.score == Double.POSITIVE_INFINITY)//node impossible, so no children
            return new ArrayList<>();
        
        ArrayList<AStarNode> ans = new ArrayList<>();
        int nextLevel = nextLevelToExpand(curNode.nodeAssignments);
        
        
        int firstIndex = nextLevel*maxNumRCs;
        for (int i=0; i<maxNumRCs; i++) {
        	int rc = unprunedRCsAtPos[firstIndex + i];
        	if (rc < 0) {
        		break;
        	}
        	
            int[] childConf = curNode.nodeAssignments.clone();
            childConf[nextLevel] = rc;
            AStarNode childNode = new AStarNode(childConf, scoreConf(childConf), useRefinement);
            ans.add(childNode);
        }
        
        return ans;
    }


    @Override
    public AStarNode rootNode() {
        //no residues assigned, so all -1's
        int[] conf = new int[numPos];
        Arrays.fill(conf,-1);
        
        AStarNode root = new AStarNode(conf, scoreConf(conf), useRefinement);
        return root;
    }
    

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        for(int rc : node.nodeAssignments){
            if(rc<0)//not fully assigned
                return false;
        }
        
        return true;
    }
    
    
    
    //operations supporting special features like dynamic A*
    
    public int nextLevelToExpand(int[] partialConf){
        //given a partially defined conformation, what level should be expanded next?
        
        if(useDynamicAStar){
            
            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;
            
            for(int level=0; level<numPos; level++){
                if(partialConf[level]<0){//position isn't already all expanded
                    
                    double levelScore = scoreExpansionLevel(level,partialConf);

                    if(levelScore>bestLevelScore){//higher score is better
                        bestLevelScore = levelScore;
                        bestLevel = level;
                    }
                }
            }
            
            if(bestLevel==-1)
                throw new RuntimeException("ERROR: No next expansion level found for dynamic A*");
            
            return bestLevel;
        }
        else {//static ordering.  
            //Let's only support the traditional ordering since dynamic will beat static for improved orderings.
            for(int level=0; level<numPos; level++){
                if(partialConf[level]<0)
                    return level;
            }
            
            throw new RuntimeException("ERROR: Can't find next expansion level for fully defined conformation");
        }
        
    }
    
    
    double scoreExpansionLevel(int level, int[] partialConf) {
        //Score expansion at the indicated level for the given partial conformation
        //for use in dynamic A*.  Higher score is better.
        
        // backup partialConf (to the stack) before we change it
        int confBak = partialConf[level];
        
        //best performing score is just 1/(sum of reciprocals of score rises for child nodes)
        double parentScore = scoreConf(partialConf);
        double reciprocalSum = 0;
        int firstIndex = level*maxNumRCs;
        for (int i=0; i<maxNumRCs; i++) {
        	int rc = unprunedRCsAtPos[firstIndex + i];
        	if (rc < 0) {
        		break;
        	}
            partialConf[level] = rc;
            double childScore = scoreConf(partialConf);
            reciprocalSum += 1.0/( childScore - parentScore );
        }
        
        // restore partialConf from the backup
        partialConf[level] = confBak;
        
        return 1.0/reciprocalSum;
    }
    
    
        
	double scoreConf(int[] partialConf) {
		
		// NOTE: might want to implement this as subclass or compose with other object
		// instead of adding a big switch here
		if(!traditionalScore) {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
		}
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some things to the stack =)
		int numPos = this.numPos;
		EnergyMatrix emat = this.emat;
		int[] unprunedRCsAtPos = this.unprunedRCsAtPos;
		int maxNumRCs = this.maxNumRCs;
		
		// TODO: could optimize memory here?
		RCTuple rcTuple = new RCTuple(partialConf);
		double score = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);//"g-score"
		
		//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
		//plus contributions associated with each of the undefined res ("h-score")
		for (int level1=0; level1<numPos; level1++) {
			
			// skip defined levels
			if (partialConf[level1] >= 0) {
				continue;
			}
			
			double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
			//resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
			
			int firstIndex = level1*maxNumRCs;
			for (int i=0; i<maxNumRCs; i++) {
				int rc1 = unprunedRCsAtPos[firstIndex + i];
				if (rc1 < 0) {
					break;
				}
				
				// OPTIMIZATION: manually inlining this is significantly slower
				// maybe it causes too much register pressure
				double rcContrib = RCContributionLB(level1, rc1, partialConf);
				
				resContribLB = Math.min(resContribLB, rcContrib);
			}
		
			score += resContribLB;
		}
		
		return score;
    }
	
	double RCContributionLB(int level1, int rc1, int[] partialConf) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some things to the stack =)
		int numPos = this.numPos;
		EnergyMatrix emat = this.emat;
		int[] unprunedRCsAtPos = this.unprunedRCsAtPos;
		int maxNumRCs = this.maxNumRCs;
		
		// OPTIMIZATION: don't even check higher terms if the energy matrix doesn't have any
		// this does wonders to CPU cache performance!
		boolean useHigherOrderTerms = emat.hasHigherOrderTerms();
		
		double rcContrib = emat.getOneBody(level1, rc1);
		
		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		for (int level2=0; level2<numPos; level2++) {
			
			int rc2 = partialConf[level2];
			boolean isDefined = rc2 >= 0;
			
			// skip this level if it's undefined and it's >= level1
			if (!isDefined && level2 >= level1) {
				continue;
			}
			
			double levelBestE = Double.POSITIVE_INFINITY;//best pairwise energy
			if (isDefined) {
				
				// consider just the one defined conformation
				double pairwiseEnergy = emat.getPairwise(level1, rc1, level2, rc2);
				if (useHigherOrderTerms) {
					//add higher-order terms that involve rc, rc2, and parts of partialConf
					//besides that only residues in definedTuple or levels below level2
					pairwiseEnergy += higherOrderContribLB(partialConf, level1, rc1, level2, rc2);
				}
				
				levelBestE = Math.min(levelBestE, pairwiseEnergy);
				
			} else {
				
				// consider all possible conformations
				int firstIndex = level2*maxNumRCs;
				for (int i=0; i<maxNumRCs; i++) {
					rc2 = unprunedRCsAtPos[firstIndex + i];
					if (rc2 < 0) {
						break;
					}
					
					double pairwiseEnergy = emat.getPairwise(level1, rc1, level2, rc2);
					if (useHigherOrderTerms) {
						pairwiseEnergy += higherOrderContribLB(partialConf, level1, rc1, level2, rc2);
					}
					
					levelBestE = Math.min(levelBestE, pairwiseEnergy);
				}
			}
			rcContrib += levelBestE;
		}
		
		return rcContrib;
	}
	
    ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf){
        //What RCs are allowed at the specified level (i.e., position num) in the given partial conf?
        ArrayList<Integer> allowedRCs = new ArrayList<>();
        if(partialConf[level]==-1) {//position undefined: consider all RCs
			int firstIndex = level*maxNumRCs;
			for (int i=0; i<maxNumRCs; i++) {
				int rc = unprunedRCsAtPos[firstIndex + i];
				if (rc < 0) {
					break;
				}
				allowedRCs.add(rc);
			}
        } else if(partialConf[level]>=0){
            allowedRCs.add(partialConf[level]);
        }
        else
            throw new UnsupportedOperationException("ERROR: Partially assigned position not yet supported in A*");
        
        return allowedRCs;
    }

    
    double higherOrderContribLB(int[] partialConf, int pos1, int rc1, int pos2, int rc2){
        //higher-order contribution for a given RC pair, when scoring a partial conf
        
        HigherTupleFinder<Double> htf = emat.getHigherOrderTerms(pos1,rc1,pos2,rc2);
        
        if(htf==null)
            return 0;//no higher-order interactions
        else
            return higherOrderContribLB(partialConf, htf, pos2);
    }
    
    
    double higherOrderContribLB(int[] partialConf, HigherTupleFinder<Double> htf, int level2){
        //recursive function to get lower bound on higher-than-pairwise terms
        //this is the contribution to the lower bound due to higher-order interactions
        //of the RC tuple corresponding to htf with "lower-numbered" residues (numbering as in scoreConf:
        //these are residues that are fully defined in partialConf, or are actually numbered <level2)

        double contrib = 0;
                
        for(int iPos : htf.getInteractingPos() ){//position has higher-order interaction with tup
            if(posComesBefore(iPos,level2,partialConf)){//interaction in right order
                //(want to avoid double-counting)
                
                double levelBestE = Double.POSITIVE_INFINITY;//best value of contribution
                //from tup-iPos interaction
                ArrayList<Integer> allowedRCs = allowedRCsAtLevel(iPos,partialConf);
                
                for( int rc : allowedRCs ){
                    
                    double interactionE = htf.getInteraction(iPos, rc);
                    
                    //see if need to go up to highers order again...
                    HigherTupleFinder<Double> htf2 = htf.getHigherInteractions(iPos, rc);
                    if(htf2!=null){
                        interactionE += higherOrderContribLB(partialConf, htf2, iPos);
                    }
                    
                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE,interactionE);
                }

                contrib += levelBestE;//add up contributions from different interacting positions iPos
            }
        }
        
        return contrib;
    }
    
    
    private boolean posComesBefore(int pos1, int pos2, int partialConf[]){
        //for purposes of contributions to traditional conf score, 
        //we go through defined and then through undefined positions (in partialConf);
        //within each of these groups we go in order of position number
        if(partialConf[pos2]>=0){//pos2 defined
            return (pos1<pos2 && partialConf[pos1]>=0);//pos1 must be defined to come before pos2
        }
        else//pos1 comes before pos2 if it's defined, or if pos1<pos2
            return (pos1<pos2 || partialConf[pos1]>=0);
    }
    
    /*
    @Override
    boolean canPruneNode(AStarNode node){
        check seq dev from wt;
    }
    
    
    
    @Override
    void refineScore(AStarNode node){//e.g. add the EPIC contribution
        node.score = betterScore();//or this could be a good place for MPLP or sthg
    }
    */
    
    
     @Override
    void refineScore(AStarNode node){
        
        if(epicMat==null)
            throw new UnsupportedOperationException("ERROR: Trying to call refinement w/o EPIC matrix");
            //later can do MPLP, etc. here
        
        if(minPartialConfs || isFullyAssigned(node))
            node.score += epicMat.minContE(node.nodeAssignments);
        
        node.scoreNeedsRefinement = false;
    }
     
     
     
    //this function computes the minimum over all full conf E's consistent with partialConf
    //for debugging only of course
    double exhaustiveScore(int[] partialConf){
        for(int pos=0; pos<partialConf.length; pos++){
            if(partialConf[pos]==-1){
                //recurse to get all options
                double score = Double.POSITIVE_INFINITY;
                for(int rc : allowedRCsAtLevel(pos,partialConf)){
                    int partialConf2[] = partialConf.clone();
                    partialConf2[pos] = rc;
                    score = Math.min(score,exhaustiveScore(partialConf2));
                }
                return score;
            }
        }
        //if we get here, conf fully defined
        return emat.getInternalEnergy( new RCTuple(partialConf) );
    }
}
