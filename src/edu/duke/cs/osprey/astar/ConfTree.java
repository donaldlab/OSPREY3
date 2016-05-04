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
    
    
    int[][] unprunedRCsAtPos;
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue
    
    
    
    
    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    boolean traditionalScore = true;
    boolean useRefinement = false;//refine nodes (might want EPIC, MPLP, or something else)
    
    boolean useDynamicAStar = true;

    
    EPICMatrix epicMat = null;//to use in refinement
    ConfSpace confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined
    
    int[] precomputedMinOffsets;
    double[] precomputedMins;
    
    // temp storage
    // NOTE: this temp storage makes this class not thread-safe!
    // but since the workload is memory-bound anyway, there isn't much benefit to parallelism
    private RCTuple rcTuple;
    int[] definedPos;
    int[] definedConf;
    int[] undefinedPos;
    
    
    public ConfTree(SearchProblem sp){
        this(sp, sp.pruneMat, sp.useEPIC);
    }
    
    public ConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC){
        numPos = sp.confSpace.numPos;
        
        // allocate temp space
        rcTuple = new RCTuple();
        definedPos = new int[numPos];
        definedConf = new int[numPos];
        undefinedPos = new int[numPos];
        
        // see which RCs are unpruned and thus available for consideration
        // pack them into an efficient int matrix
        unprunedRCsAtPos = new int[numPos][];
        for (int pos=0; pos<numPos; pos++) {
        	ArrayList<Integer> srcRCs = pruneMat.unprunedRCsAtPos(pos);
        	int[] destRCs = new int[srcRCs.size()];
        	for (int i=0; i<srcRCs.size(); i++) {
        		destRCs[i] = srcRCs.get(i);
        	}
        	unprunedRCsAtPos[pos] = destRCs;
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
        
        // OPTIMIZATION: precompute some mins to speed up scoreConf()
        // only works if we're not using higher order terms though
        precomputedMinOffsets = null;
        precomputedMins = null;
        if (!emat.hasHigherOrderTerms()) {
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
    }
    
    
    public BigInteger getNumConformations() {
    	BigInteger num = BigInteger.valueOf(1);
    	for (int pos=0; pos<numPos; pos++) {
    		num = num.multiply(BigInteger.valueOf(unprunedRCsAtPos[pos].length));
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
        
        
        for (int rc : unprunedRCsAtPos[nextLevel]) {
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
        for (int rc : unprunedRCsAtPos[level]) {
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
		
		rcTuple.set(partialConf);
		double gscore = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);
		
		// OPTIMIZATION: separating scoring by defined vs undefined residues noticeably improves CPU cache performance
		// (compared to scoring in order of residue and alternating between the two different types of scoring)
		
		int[] definedPos = this.definedPos;
		int[] definedConf = this.definedConf;
		int[] undefinedPos = this.undefinedPos;
		
		// split conformation into defined and undefined residues
		int numDefined = 0;
		for (int pos=0; pos<partialConf.length; pos++) {
			if (partialConf[pos] >= 0) {
				numDefined++;
			}
		}
		int numUndefined = partialConf.length - numDefined;
		int definedIndex = 0;
		int undefinedIndex = 0;
		for (int pos=0; pos<partialConf.length; pos++) {
			int rc = partialConf[pos];
			if (rc >= 0) {
				definedPos[definedIndex] = pos;
				definedConf[definedIndex] = rc;
				definedIndex++;
			} else {
				undefinedPos[undefinedIndex] = pos;
				undefinedIndex++;
			}
		}
		
		return gscore + getHScore(partialConf, numDefined, numUndefined);
    }
	
	private double getHScore(int[] conf, int numDefined, int numUndefined) {
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some things to the stack =)
		int[][] unprunedRCsAtPos = this.unprunedRCsAtPos;
		int[] undefinedPos = this.undefinedPos;
		
		//score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
		//plus contributions associated with each of the undefined res ("h-score")
		double hscore = 0;
		for (int i=0; i<numUndefined; i++) {
			int pos1 = undefinedPos[i];
			
			//lower bound on contribution of this residue
			//resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
			double resContribLB = Double.POSITIVE_INFINITY;
			for (int rc1 : unprunedRCsAtPos[pos1]) {
				// OPTIMIZATION: manually inlining this is noticeably slower
				// maybe it causes too much register pressure
				double rcContrib = RCContributionLB(conf, numDefined, numUndefined, pos1, rc1);
				resContribLB = Math.min(resContribLB, rcContrib);
			}
		
			hscore += resContribLB;
		}
		return hscore;
	}
	
	double RCContributionLB(int[] conf, int numDefined, int numUndefined, int pos1, int rc1) {
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some references to the stack =)
		EnergyMatrix emat = this.emat;
		int[][] unprunedRCsAtPos = this.unprunedRCsAtPos;
		int[] precomputedMinOffsets = this.precomputedMinOffsets;
		double[] precomputedMins = this.precomputedMins;
		int[] definedPos = this.definedPos;
		int[] definedConf = this.definedConf;
		int[] undefinedPos = this.undefinedPos;
		
		// OPTIMIZATION: don't even check higher terms if the energy matrix doesn't have any
		// this does wonders to CPU cache performance!
		boolean useHigherOrderTerms = emat.hasHigherOrderTerms();
		
		double rcContrib = emat.getOneBody(pos1, rc1);
		
		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		
		// first pass, defined residues
		for (int i=0; i<numDefined; i++) {
			int pos2 = definedPos[i];
			int rc2 = definedConf[i];
			
			rcContrib += emat.getPairwise(pos1, rc1, pos2, rc2);
			if (useHigherOrderTerms) {
				//add higher-order terms that involve rc, rc2, and parts of partialConf
				//besides that only residues in definedTuple or levels below pos2
				rcContrib += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
			}
		}
		
		// second pass, undefined residues
		for (int i=0; i<numUndefined; i++) {
			int pos2 = undefinedPos[i];
			if (pos2 >= pos1) {
				break;
			}
				
			// if we're using higher order terms, we have to compute the mins here
			if (useHigherOrderTerms) {
				
				// min over all possible conformations
				double minEnergy = Double.POSITIVE_INFINITY;
				for (int rc2 : unprunedRCsAtPos[pos2]) {
					double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
					pairwiseEnergy += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
					minEnergy = Math.min(minEnergy, pairwiseEnergy);
				}
				rcContrib += minEnergy;
				
			// but otherwise, we can use the precomputations
			} else {
				int index = precomputedMinOffsets[pos1*(pos1 - 1)/2 + pos2] + rc1;
				rcContrib += precomputedMins[index];
			}
		}
		
		return rcContrib;
	}
	
    ArrayList<Integer> allowedRCsAtLevel(int level, int[] partialConf){
        //What RCs are allowed at the specified level (i.e., position num) in the given partial conf?
        ArrayList<Integer> allowedRCs = new ArrayList<>();
        if(partialConf[level]==-1) {//position undefined: consider all RCs
        	for (int rc : unprunedRCsAtPos[level]) {
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
