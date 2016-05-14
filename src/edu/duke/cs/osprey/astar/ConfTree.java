/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

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

    protected int numPos;
    protected EnergyMatrix emat;
    
    
    protected int[][] unprunedRCsAtPos;
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue
    
    
    
    
    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    protected boolean traditionalScore = true;
    protected boolean useRefinement = false;//refine nodes (might want EPIC, MPLP, or something else)
    
    protected boolean useDynamicAStar = true;

    
    EPICMatrix epicMat = null;//to use in refinement
    protected ConfSpace confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined
    
    // temp storage
    // NOTE: this temp storage makes this class not thread-safe!
    // but since the workload is memory-bound anyway, there isn't much benefit to parallelism
    protected RCTuple rcTuple;
    protected int numDefined;
    protected int numUndefined;
    protected int[] definedPos;
    protected int[] definedRCs;
    protected int[] undefinedPos;
    private int[] childConf;
    
    
    public ConfTree(SearchProblem sp){
        this(sp, sp.pruneMat, sp.useEPIC);
    }
    
    public ConfTree(SearchProblem sp, PruningMatrix pruneMat, boolean useEPIC){
    	
		// NOTE: might want to implement this as subclass or compose with other object
		// instead of adding a big switch here
		if(!traditionalScore) {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
		}
		
        numPos = sp.confSpace.numPos;
        
        // allocate temp space
        rcTuple = new RCTuple();
        numDefined = 0;
        numUndefined = 0;
        definedPos = new int[numPos];
        definedRCs = new int[numPos];
        undefinedPos = new int[numPos];
        childConf = new int[numPos];
        
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
    }
    
    
    @Override
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
        
        if(curNode.getScore() == Double.POSITIVE_INFINITY)//node impossible, so no children
            return new ArrayList<>();
        
        ArrayList<AStarNode> ans = new ArrayList<>();
        int nextLevel = nextLevelToExpand(curNode);
        
        
        for (int rc : unprunedRCsAtPos[nextLevel]) {
            int[] childConf = curNode.getNodeAssignments().clone();
            childConf[nextLevel] = rc;
            AStarNode childNode = new AStarNode(childConf, useRefinement);
            scoreNodeDifferential(curNode, childNode, nextLevel, rc);
            ans.add(childNode);
        }
        
        return ans;
    }


	@Override
	public AStarNode rootNode() {
		
		//no residues assigned, so all -1's
		int[] conf = new int[numPos];
		Arrays.fill(conf, -1);
		
		AStarNode root = new AStarNode(conf, useRefinement);
		scoreNode(root);
		return root;
	}
    

    @Override
    public boolean isFullyAssigned(AStarNode node) {
        for(int rc : node.getNodeAssignments()){
            if(rc<0)//not fully assigned
                return false;
        }
        
        return true;
    }
    
    
    
    //operations supporting special features like dynamic A*
    
    public int nextLevelToExpand(AStarNode parentNode) {
        //given a partially defined conformation, what level should be expanded next?
        
        int[] conf = parentNode.getNodeAssignments();
        if(useDynamicAStar){
            
            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;
            
            for(int level=0; level<numPos; level++){
                if(conf[level]<0){//position isn't already all expanded
                    
                    double levelScore = scoreExpansionLevel(parentNode, level);

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
                if(conf[level]<0)
                    return level;
            }
            
            throw new RuntimeException("ERROR: Can't find next expansion level for fully defined conformation");
        }
        
    }
    
    
    double scoreExpansionLevel(AStarNode parentNode, int level) {
        //Score expansion at the indicated level for the given partial conformation
        //for use in dynamic A*.  Higher score is better.
    	
        //best performing score is just 1/(sum of reciprocals of score rises for child nodes)
        double parentScore = parentNode.getScore();
        double reciprocalSum = 0;
        for (int rc : unprunedRCsAtPos[level]) {
            double childScore = scoreConfDifferential(parentNode, level, rc);
            reciprocalSum += 1.0/( childScore - parentScore );
        }
        
        return 1.0/reciprocalSum;
    }
    
    // NOTE about scoring:
	// score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
	// plus contributions associated with each of the undefined res ("h-score")
    
    // we'll optimize that calculation by treating defined and undefined positions separately though
    // this improves CPU cache performance, which is important since the big energy matrix makes
    // this computation mostly memory-bound
    
    protected void splitPositions(int[] conf) {
    	
    	int numPos = this.numPos;
		int[] definedPos = this.definedPos;
		int[] definedRCs = this.definedRCs;
		int[] undefinedPos = this.undefinedPos;
		
		// split conformation into defined and undefined residues
		int numDefined = 0;
		for (int pos=0; pos<numPos; pos++) {
			if (conf[pos] >= 0) {
				numDefined++;
			}
		}
		int definedIndex = 0;
		int undefinedIndex = 0;
		for (int pos=0; pos<numPos; pos++) {
			int rc = conf[pos];
			if (rc >= 0) {
				definedPos[definedIndex] = pos;
				definedRCs[definedIndex] = rc;
				definedIndex++;
			} else {
				undefinedPos[undefinedIndex] = pos;
				undefinedIndex++;
			}
		}
		
		this.numDefined = numDefined;
		this.numUndefined = numPos - numDefined;
		
		this.rcTuple.set(conf);
    }
    
    protected void resetSplitPositions() {
    	numDefined = 0;
    	numUndefined = 0;
    }
    
    protected void assertSplitPositions() {
    	assert (numDefined + numUndefined == numPos) :
    		"call splitPostions(conf) before calling this function!";
    }
    
    protected double scoreNodeDifferential(AStarNode parent, AStarNode child, int nextPos, int nextRc) {
    	return scoreNode(child);
    }
    
    protected double scoreNode(AStarNode node) {
    	
    	// just route to scoreConf(), subclasses can do fancier things
    	double score = scoreConf(node.getNodeAssignments());
    	node.setScore(score);
    	return score;
    }
    
    protected double scoreConfDifferential(AStarNode parentNode, int nextPos, int nextRc) {
    	
    	// just route to scoreConf(), subclasses can do fancier things
    	System.arraycopy(parentNode.getNodeAssignments(), 0, childConf, 0, numPos);
    	assert (childConf[nextPos] < 0);
       	childConf[nextPos] = nextRc;
       	return scoreConf(childConf);
    }
    
	protected double scoreConf(int[] conf) {
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some things to the stack =)
		int[][] unprunedRCsAtPos = this.unprunedRCsAtPos;
		int[] undefinedPos = this.undefinedPos;
		
		splitPositions(conf);
		
		// compute g-score
		rcTuple.set(conf);
    	double gscore = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);
		
		// compute h-score
		double hscore = 0;
		for (int i=0; i<numUndefined; i++) {
			int pos1 = undefinedPos[i];
			
			//lower bound on contribution of this residue
			//resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
			double resContribLB = Double.POSITIVE_INFINITY;
			int[] rc1s = unprunedRCsAtPos[pos1];
			int n1 = rc1s.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rc1s[j];
				// OPTIMIZATION: manually inlining this is noticeably slower
				// maybe it causes too much register pressure
				double rcContrib = getUndefinedRCEnergy(conf, pos1, rc1, j);
				resContribLB = Math.min(resContribLB, rcContrib);
			}
		
			hscore += resContribLB;
		}
		
		resetSplitPositions();
		
		return gscore + hscore;
    }
	
	double getUndefinedRCEnergy(int[] conf, int pos1, int rc1, int rc1i) {
		assertSplitPositions();
		//Provide a lower bound on what the given rc at the given level can contribute to the energy
		//assume partialConf and definedTuple
		
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations (like preferring stack over heap) can make an impact
		
		// that said, let's copy some references to the stack =)
		EnergyMatrix emat = this.emat;
		int numDefined = this.numDefined;
		int numUndefined = this.numUndefined;
		int[] definedPos = this.definedPos;
		int[] definedRCs = this.definedRCs;
		int[] undefinedPos = this.undefinedPos;
		
		double rcContrib = emat.getOneBody(pos1, rc1);
		
		//for this kind of lower bound, we need to split up the energy into the defined-tuple energy
		//plus "contributions" for each undefined residue
		//so we'll say the "contribution" consists of any interactions that include that residue
		//but do not include higher-numbered undefined residues
		
		// first pass, defined residues
		for (int i=0; i<numDefined; i++) {
			int pos2 = definedPos[i];
			int rc2 = definedRCs[i];
			
			rcContrib += emat.getPairwise(pos1, rc1, pos2, rc2);
			//add higher-order terms that involve rc, rc2, and parts of partialConf
			//besides that only residues in definedTuple or levels below pos2
			rcContrib += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
		}
		
		// second pass, undefined residues
		for (int i=0; i<numUndefined; i++) {
			int pos2 = undefinedPos[i];
			if (pos2 >= pos1) {
				break;
			}
			
			rcContrib += getMinPairwiseEnergy(conf, pos1, rc1, rc1i, pos2);
		}
		
		return rcContrib;
	}
	
	protected double getMinPairwiseEnergy(int[] conf, int pos1, int rc1, int rc1i, int pos2) {
		
		EnergyMatrix emat = this.emat;
		
		// min over all possible conformations
		double minEnergy = Double.POSITIVE_INFINITY;
		for (int rc2 : this.unprunedRCsAtPos[pos2]) {
			double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
			pairwiseEnergy += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
			minEnergy = Math.min(minEnergy, pairwiseEnergy);
		}
		return minEnergy;
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
            node.setScore(node.getScore() + epicMat.minContE(node.getNodeAssignments()));
        
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
