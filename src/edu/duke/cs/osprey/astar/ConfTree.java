/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.confspace.SearchSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
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
    SearchSpace searchSpace;
    
    
    
    ArrayList<ArrayList<Integer>> unprunedRCsAtPos = new ArrayList<>();
    //get from searchSpace when initializing!
    //These are lists of residue-specific RC numbers for the unpruned RCs at each residue
    
    
    
    
    //ADVANCED SCORING METHODS: TO CHANGE LATER (EPIC, MPLP, etc.)
    boolean traditionalScore = true;
    boolean useRefinement = false;//refine nodes (might want EPIC, MPLP, or something else)
    
    boolean useDynamicAStar = false;

    
    
    public ConfTree(SearchSpace acs) {
        searchSpace = acs;
        numPos = searchSpace.confSpace.numPos;
        
        for(int pos=0; pos<numPos; pos++){
            unprunedRCsAtPos.add( searchSpace.pruneMat.unprunedRCsAtPos(pos) );
        }
    }
    
    
    
    
    @Override
    ArrayList<AStarNode> getChildren(AStarNode curNode) {
        
        if(isFullyAssigned(curNode))
            throw new RuntimeException("ERROR: Can't expand a fully assigned A* node");
        
        ArrayList<AStarNode> ans = new ArrayList<>();
        int nextLevel = nextLevelToExpand(curNode.nodeAssignments);
        
        for(int rc : unprunedRCsAtPos.get(nextLevel) ){
            int[] childConf = curNode.nodeAssignments.clone();
            childConf[nextLevel] = rc;
            AStarNode childNode = new AStarNode(childConf, scoreConf(childConf), useRefinement);
            ans.add(childNode);
        }
        
        return ans;
    }


    @Override
    AStarNode rootNode() {
        //no residues assigned, so all -1's
        int[] conf = new int[numPos];
        Arrays.fill(conf,-1);
        
        AStarNode root = new AStarNode(conf, scoreConf(conf), useRefinement);
        return root;
    }
    

    @Override
    boolean isFullyAssigned(AStarNode node) {
        for(int rc : node.nodeAssignments){
            if(rc<0)//not fully assigned
                return false;
        }
        
        return true;
    }
    
    
    
    //operations supporting special features like dynamic A*
    
    int nextLevelToExpand(int[] partialConf){
        //given a partially defined conformation, what level should be expanded next?
        
        if(useDynamicAStar){
            
            int bestLevel = -1;
            double bestLevelScore = Double.NEGATIVE_INFINITY;
            
            for(int level=0; level<numPos; level++){
                double levelScore = scoreExpansionLevel(level,partialConf);
                
                if(levelScore>bestLevelScore){//higher score is better
                    bestLevelScore = levelScore;
                    bestLevel = level;
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
    
    
    double scoreExpansionLevel(int level, int[] partialConf){
        //Score expansion at the indicated level for the given partial conformation
        //for use in dynamic A*.  Higher score is better.
        
        //best performing score is just 1/(sum of reciprocals of score rises for child nodes)
        double parentScore = scoreConf(partialConf);
        int[] expandedConf = partialConf.clone();
        
        double reciprocalSum = 0;
        
        for(int rc : unprunedRCsAtPos.get(level) ){
            expandedConf[level] = rc;
            double childScore = scoreConf(expandedConf);
            
            reciprocalSum += 1.0 / ( childScore - parentScore );
        }
        
        double score = 1. / reciprocalSum;
        
        return score;
    }
    
    
    
    double scoreConf(int[] partialConf){
        if(traditionalScore){
            RCTuple definedTuple = new RCTuple(partialConf);
            
            double score = searchSpace.emat.getInternalEnergy( definedTuple );//"g-score"
            
            //score works by breaking up the full energy into the energy of the defined set of residues ("g-score"),
            //plus contributions associated with each of the undefined res ("h-score")
            for(int level=0; level<numPos; level++){
                if(partialConf[level]<0){//level not fully defined
                    
                    double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
                    //resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
                    
                    for ( int rc : unprunedRCsAtPos.get(level) ) {
                        resContribLB = Math.min(resContribLB, RCContributionLB(level,rc,definedTuple,partialConf));
                    }
                
                    score += resContribLB;
                }
            }
            
            return score;
        }
        else {
            //other possibilities include MPLP, etc.
            //But I think these are better used as refinements
            //we may even want multiple-level refinement
            throw new RuntimeException("Advanced A* scoring methods not implemented yet!");
        }
    }
    
    
    
    double RCContributionLB(int level, int rc, RCTuple definedTuple, int[] partialConf){
        //Provide a lower bound on what the given rc at the given level can contribute to the energy
        //assume partialConf and definedTuple
        
        double rcContrib = searchSpace.emat.getOneBody(level,rc);
        
        //for this kind of lower bound, we need to split up the energy into the defined-tuple energy
        //plus "contributions" for each undefined residue
        //so we'll say the "contribution" consists of any interactions that include that residue
        //but do not include higher-numbered undefined residues
        for(int level2=0; level2<level; level2++){
            
            if(partialConf[level2]>=0 || level2<level){//lower-numbered or defined residues
                
                double levelBestE = Double.POSITIVE_INFINITY;//best pairwise energy
                
                ArrayList<Integer> allowedRCs;
                if(partialConf[level2]==-1)//position undefined: consider all RCs
                    allowedRCs = unprunedRCsAtPos.get(level2);
                else if(partialConf[level2]>=0){
                    allowedRCs = new ArrayList<>();
                    allowedRCs.add(partialConf[level2]);
                }
                else
                    throw new UnsupportedOperationException("ERROR: Partially assigned position not yet supported in A*");

                
                for( int rc2 : allowedRCs ){
                    
                    double interactionE = searchSpace.emat.getPairwise(level,rc,level2,rc2);
                    
                    //double higherLB = higherOrderContribLB(partialConf,level,rc,level2,rc2,);
                    //add higher-order terms that involve rc, rc2, and
                    
                    //interactionE += higherLB;
                    
                    //besides that only residues in definedTuple or levels below level2
                    levelBestE = Math.min(levelBestE,interactionE);
                }

                rcContrib += levelBestE;
            }
        }
        
        return rcContrib;
    }
    
    /*
    double higherOrderContribLB(){
        //recursive function to get lower bound, in min-sum-min fashion, on triples+ terms
    }
    
    
    
    @Override
    boolean canPruneNode(AStarNode node){
        check seq dev from wt;
    }
    
    
    
    
    @Override
    void refineScore(AStarNode node){//e.g. add the EPIC contribution
        node.score = betterScore();//or this could be a good place for MPLP or sthg
    }
    */
}
