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
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.NewEPICMatrix;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * @author mhall44
 */
public class ConfTree<T extends AStarNode> extends AStarTree<T> {
    //This implementation of an A* tree is intended for conformational search
    //AStarNode.nextAssignment is an array of length numPos; each position
    //stores the assigned RC, or -1 to indicate an unassigned position
    //this class supports both "traditional" (static, simple heuristic) A*
    //and improvements like dynamic A*
    //we may also want to allow other negative indices, to indicate partially assigned RCs
	
	private static final long serialVersionUID = -2776047349227703884L;


	public static ConfTree<FullAStarNode> makeFull(SearchProblem search) {
		return makeFull(search, search.pruneMat);
	}
	
	public static ConfTree<FullAStarNode> makeFull(SearchProblem search, PruningMatrix pmat) {
		return new ConfTree<FullAStarNode>(new FullAStarNode.Factory(search.confSpace.numPos), search, pmat);
	}
        
        public static ConfTree<FullAStarNode> makeFull(SearchProblem search, PruningMatrix pmat, GMECMutSpace gms) {
		return new ConfTree<FullAStarNode>(new FullAStarNode.Factory(search.confSpace.numPos), search, pmat, search.useEPIC, gms);
	}
        
        public static ConfTree<FullAStarNode> makeFull(PrecomputedMatrices precompMat, GMECMutSpace gms,
                boolean useTupExpForSearch, boolean useEPIC, EPICSettings epicSettings, int numPos){
            return new ConfTree<FullAStarNode>(new FullAStarNode.Factory(numPos),
                    useTupExpForSearch, precompMat.getEmat(), null, precompMat.getEpicMat(), precompMat.getLuteMat(),
                    precompMat.getPruneMat(), useEPIC, epicSettings, gms);
        }
	
	private AStarNode.Factory<T> nodeFactory;

    protected int numPos;
    protected EnergyMatrix emat;
    
    
    protected int[][] unprunedRCsAtPos;
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue
    
    
    
    
    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    protected boolean traditionalScore = true;
    protected boolean useRefinement = false;//refine nodes (might want EPIC, MPLP, or something else)
    
    protected boolean useDynamicAStar = true;

    //We can have at most one of these
    protected EPICMatrix epicMat = null;//to use in refinement
    protected NewEPICMatrix newEPICMat = null;
    //protected ConfSpace confSpace = null;//conf space to use with epicMat if we're doing EPIC minimization w/ SAPE
    protected boolean minPartialConfs = false;//whether to minimize partially defined confs with EPIC, or just fully defined
    
    // temp storage
    // NOTE: this temp storage makes this class not thread-safe!
    // but since the workload is memory-bound anyway, there isn't much benefit to parallelism
    protected RCTuple rcTuple;
    protected int numDefined;
    protected int numUndefined;
    protected int[] definedPos;
    protected int[] definedRCs;
    protected int[] undefinedPos;
    protected int[] childConf;
    
    private GMECMutSpace mutSpace;
    
    
    public ConfTree(AStarNode.Factory<T> nodeFactory, SearchProblem sp){
        this(nodeFactory, sp, sp.pruneMat, sp.useEPIC, null);
    }
    
    public ConfTree(AStarNode.Factory<T> nodeFactory, SearchProblem sp, PruningMatrix pruneMat){
        this(nodeFactory, sp, pruneMat, sp.useEPIC, null);
    }
    
    public ConfTree(AStarNode.Factory<T> nodeFactory, SearchProblem sp, PruningMatrix pruneMat, 
            boolean useEPIC, GMECMutSpace gms){
        this(nodeFactory, sp.useTupExpForSearch, sp.emat, sp.epicMat, null,  
                sp.tupExpEMat, pruneMat, sp.useEPIC, sp.epicSettings, null);
    }
    
    
    //matrices/settings we won't use can be null
    public ConfTree(AStarNode.Factory<T> nodeFactory, boolean useTupExpForSearch,
            EnergyMatrix emat, EPICMatrix epicMat, NewEPICMatrix newEPICMat, EnergyMatrix tupExpEMat, PruningMatrix pruneMat, 
            boolean useEPIC, EPICSettings epicSettings, GMECMutSpace gms){
    	
		// NOTE: might want to implement this as subclass or compose with other object
		// instead of adding a big switch here
		if(!traditionalScore) {
			//other possibilities include MPLP, etc.
			//But I think these are better used as refinements
			//we may even want multiple-level refinement
			throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
		}
		
		this.nodeFactory = nodeFactory;
		// AAO 2016: for fully defined sequences, this change has no effect. for partially defined
		// sequences, however, it obviates the need for a separate function to initialize the data
		// structures below.
        numPos = pruneMat.getNumPos();
        
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
        if(useTupExpForSearch)
            this.emat = tupExpEMat;
        else {
            this.emat = emat;
            
            if(useEPIC){//include EPIC in the search
                useRefinement = true;
                this.epicMat = epicMat;
                this.newEPICMat = newEPICMat;
                if((epicMat==null)==(newEPICMat==null)){
                    throw new RuntimeException("ERROR: to use EPIC in A* need exactly one of old and new EPIC matrices");
                }
                //confSpace = sp.confSpace;
                minPartialConfs = epicSettings.minPartialConfs;
            }
        }
        
        this.mutSpace = gms;
        if(gms!=null)//GMECMutSpace requires static ordering
            useDynamicAStar = false;
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
    public ArrayList<T> getChildren(T curNode) {
        
        if(isFullyAssigned(curNode))
            throw new RuntimeException("ERROR: Can't expand a fully assigned A* node");
        
        if(curNode.getScore() == Double.POSITIVE_INFINITY)//node impossible, so no children
            return new ArrayList<>();
        
        ArrayList<T> ans = new ArrayList<>();
        int nextLevel = nextLevelToExpand(curNode);
        
        splitPositions(curNode);
        
        for (int rc : unprunedRCsAtPos[nextLevel]) {
            if(mutSpace!=null){
                if(!mutSpace.isNewRCAllowed(curNode.getNodeAssignments(),curNode.getLevel(),rc)){
                    continue;
                }
            }
            
            T childNode = nodeFactory.make(curNode, nextLevel, rc);
            childNode.setScoreNeedsRefinement(useRefinement);
            scoreNodeDifferential(curNode, childNode, nextLevel, rc);
            ans.add(childNode);
        }
        
        resetSplitPositions();
        
        return ans;
    }
    
	@Override
	public T rootNode() {
		
		//no residues assigned, so all -1's
		int[] conf = new int[numPos];
		Arrays.fill(conf, -1);
		
		T root = nodeFactory.makeRoot();
		root.setScoreNeedsRefinement(useRefinement);
		splitPositions(root);
		scoreNode(root);
		resetSplitPositions();
		return root;
	}
    

    @Override
    public boolean isFullyAssigned(T node) {
    	return node.isFullyDefined();
    }
    
    
    
    //operations supporting special features like dynamic A*
    
    public int nextLevelToExpand(T parentNode) {
        //given a partially defined conformation, what level should be expanded next?
    	
        int bestLevel = -1;
        
    	splitPositions(parentNode);
    	
        if (useDynamicAStar) {
            
            double bestLevelScore = Double.NEGATIVE_INFINITY;
            
            for (int i=0; i<numUndefined; i++) {
            	int level = undefinedPos[i];
                    
				double levelScore = scoreExpansionLevel(parentNode, level);

				if(levelScore>bestLevelScore){//higher score is better
					bestLevelScore = levelScore;
					bestLevel = level;
                }
            }
            
        } else {//static ordering.  
            //Let's only support the traditional ordering since dynamic will beat static for improved orderings.
        	if (numUndefined > 0) {
        		bestLevel = undefinedPos[0];
        	}
        }
        	
        resetSplitPositions();
        
        if(bestLevel==-1) {
        	throw new RuntimeException("ERROR: No next expansion level found for dynamic A*");
        } else {
        	return bestLevel;
        }
    }
    
    
    double scoreExpansionLevel(T parentNode, int level) {
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
    
    protected void splitPositions(T node) {
    	
    	// make sure we're not split already
    	assert (numDefined == 0 && numUndefined == 0);
    	
    	int[] conf = node.getNodeAssignments();
    	
		// split conformation into defined and undefined residues
		numDefined = 0;
		numUndefined = 0;
		for (int pos=0; pos<numPos; pos++) {
			int rc = conf[pos];
			if (rc >= 0) {
				definedPos[numDefined] = pos;
				definedRCs[numDefined] = rc;
				numDefined++;
			} else {
				undefinedPos[numUndefined] = pos;
				numUndefined++;
			}
		}
		
		assert (numDefined + numUndefined == numPos);
    }
    
    private void resetSplitPositions() {
    	numDefined = 0;
    	numUndefined = 0;
    }
    
    protected void assertSplitPositions() {
    	assert (numDefined + numUndefined == numPos) :
    		"call splitPostions(node) before calling this function!";
    }
    
    protected double scoreNodeDifferential(T parent, T child, int childPos, int childRc) {
    	assertSplitPositions();
    	
    	// just route to scoreConfDifferential()
    	double score = scoreConfDifferential(parent, childPos, childRc);
    	child.setScore(score);
    	return score;
    }
    
    protected double scoreNode(T node) {
    	assertSplitPositions();
    	
    	// just route to scoreConfDifferential()
    	double score = scoreConfDifferential(node, -1, -1);
    	node.setScore(score);
    	return score;
    }
    
    protected double scoreConfDifferential(T parentNode, int childPos, int childRc) {
    	assertSplitPositions();
    	
		// OPTIMIZATION: this function gets hit a LOT!
		// so even really pedantic optimizations can make an impact
		
    	// get the full conf, start with the parent first
    	int[] conf = parentNode.getNodeAssignments();
    	
    	// but if this is actually a child node, switch to the child conf
    	if (childPos >= 0) {
    		
    		// parent shouldn't be assigned here
    		assert (conf[childPos] < 0);
    		
    		// make the child conf
    		System.arraycopy(conf, 0, childConf, 0, numPos);
    		childConf[childPos] = childRc;
    		conf = childConf;
    	}
		
		// compute g-score
       	rcTuple.set(conf);
    	double gscore = emat.getConstTerm() + emat.getInternalEnergy(rcTuple);
		
		// compute h-score
		double hscore = 0;
		for (int i=0; i<numUndefined; i++) {
			int pos1 = undefinedPos[i];
			
			// skip if defined in child
			if (pos1 == childPos) {
				continue;
			}
			
			//lower bound on contribution of this residue
			//resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
			double resContribLB = Double.POSITIVE_INFINITY;
			int[] rc1s = unprunedRCsAtPos[pos1];
			int n1 = rc1s.length;
			for (int j=0; j<n1; j++) {
				int rc1 = rc1s[j];
				// OPTIMIZATION: manually inlining this is noticeably slower
				// maybe it causes too much register pressure
				double rcContrib = getUndefinedRCEnergy(conf, pos1, rc1, j, childPos, childRc);
				resContribLB = Math.min(resContribLB, rcContrib);
			}
		
			hscore += resContribLB;
		}
		
		return gscore + hscore;
    }
	
	private double getUndefinedRCEnergy(int[] conf, int pos1, int rc1, int rc1i, int childPos, int childRc) {
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
			
			assert (pos2 != childPos);
			
			rcContrib += emat.getPairwise(pos1, rc1, pos2, rc2);
			//add higher-order terms that involve rc, rc2, and parts of partialConf
			//besides that only residues in definedTuple or levels below pos2
			rcContrib += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
		}
		
		// if the child has a new definition, add that too
		if (childPos >= 0) {
			rcContrib += emat.getPairwise(pos1, rc1, childPos, childRc);
			rcContrib += higherOrderContribLB(conf, pos1, rc1, childPos, childRc);
		}
		
		// second pass, undefined residues
		for (int i=0; i<numUndefined; i++) {
			int pos2 = undefinedPos[i];
			if (pos2 >= pos1) {
				break;
			}
			
			// skip if defined in child
			if (pos2 == childPos) {
				continue;
			}
			
			// min over all possible conformations
			double minEnergy = Double.POSITIVE_INFINITY;
			for (int rc2 : this.unprunedRCsAtPos[pos2]) {
				double pairwiseEnergy = emat.getPairwise(pos1, rc1, pos2, rc2);
				pairwiseEnergy += higherOrderContribLB(conf, pos1, rc1, pos2, rc2);
				minEnergy = Math.min(minEnergy, pairwiseEnergy);
			}
			
			rcContrib += minEnergy;
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
    boolean canPruneNode(T node){
        check seq dev from wt;
    }
    
    
    
    @Override
    void refineScore(T node){//e.g. add the EPIC contribution
        node.score = betterScore();//or this could be a good place for MPLP or sthg
    }
    */
    
    
     @Override
    public void refineScore(T node){
        
            
        if(minPartialConfs || isFullyAssigned(node)){
            if(epicMat==null)
                node.setScore(node.getScore() + newEPICMat.minContE(node.getNodeAssignments()));
            else {
                node.setScore(epicMat.minimizeEnergy(new RCTuple(node.getNodeAssignments()), true));
            }
        
            node.setScoreNeedsRefinement(false);
        }
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
