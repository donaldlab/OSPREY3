/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMethod.CheckSumType;
import java.util.ArrayList;

/**
 *
 * @author mhall44
 */
public class Pruner {
    //This representation based on "candidate" and "competitor" objects and iterators
    //will let us cut down a lot on the code redundancy
    //and will make it way easier to update things when we need to
    //"DEECandidate", etc. can represent either a rotamer or a pair (or a triple?)
    //the candidate, competitor, and witness iterators can be reused between most types of Pruner
    //(usually we just need an iterator over rotamers, or an iterator over pairs)
    //and then each type of Pruner has its own pruning condition
    
    SearchProblem searchSpace;
    
    boolean typeDep;//use type-dependent pruning
    double boundsThreshold;//any conformations above this threshold are liable to Bounds pruning
    
    double pruningInterval;//any conformations more than Ew above optimal are liable to competitive pruning
    //(this includes I in iMinDEE, since in this case we are using the pairwise-minimal energies
    //and "optimal" refers to the lowest pairwise bound)
    
    
    public Pruner(SearchProblem searchSpace, boolean typeDep, double boundsThreshold, double pruningInterval) {
        this.searchSpace = searchSpace;
        this.typeDep = typeDep;
        this.boundsThreshold = boundsThreshold;
        this.pruningInterval = pruningInterval;
    }
    
    
    boolean prune(String methodName){
        //convenience method
        return prune(PruningMethod.getMethod(methodName));
    }
    
    
    boolean prune(PruningMethod method){
        //return whether we pruned something or not
        //A PruningMethod will just contain settings specifying the kind of Pruner
        //and return the right kinds of iterators and evaluate pruning conditions
        
        System.out.println("Starting pruning with " + method.name());
        
        boolean prunedSomething = false;
        boolean prunedSomethingThisCycle;

        do {
            prunedSomethingThisCycle = false;
            
            ArrayList<RCTuple> candidates = enumerateAllCandidates(method.numPos);
            //we can have classes RotIterator and RotPairIterator

            for( RCTuple cand : candidates ){
                
                if( ! searchSpace.pruneMat.isPruned(cand) ){
                    
                    boolean prunedCandidate = false;
                    
                    if(method.useCompetitor){
                        //any of the candidates may also be used as a competitor...
                        for(RCTuple competitor : searchSpace.pruneMat.unprunedRCTuplesAtPos(cand.pos)){
                            
                            if(cand.isSameTuple(competitor))//you can't pruned yourself
                                continue;
                            
                            if(typeDep){
                                //can only prune using competitors of same res type
                                if( !resTypesMatch(cand,competitor) )
                                    continue;
                            }
                                                        
                            if( canPrune(cand,competitor,method.cst) ){
                                prunedCandidate = true;
                                break;
                            }
                        }
                    }
                    else {
                        prunedCandidate = canPrune(cand,method.cst);//non-competitive pruning attempt
                    }
                    
                    if(prunedCandidate){
                        searchSpace.pruneMat.markAsPruned(cand);
                        prunedSomething = true;
                        prunedSomethingThisCycle = true;
                    }
                }
            }
            
        } while(prunedSomethingThisCycle);
        
        return prunedSomething;
    }
    
    
    
    boolean resTypesMatch(RCTuple tup1, RCTuple tup2){
        //Do the tuples have the same residue types?  If they do then we can
        //use one to prune the other by type-dependent DEE
        int numPosInTup = tup1.pos.size();
        
        for(int indexInTup=0; indexInTup<numPosInTup; indexInTup++){
            
            int pos1 = tup1.pos.get(indexInTup);
            int rc1 = tup1.RCs.get(indexInTup);
            String type1 = searchSpace.confSpace.posFlex.get(pos1).RCs.get(rc1).AAType;
            
            int pos2 = tup2.pos.get(indexInTup);
            int rc2 = tup2.RCs.get(indexInTup);
            String type2 = searchSpace.confSpace.posFlex.get(pos2).RCs.get(rc2).AAType;
            
            if(!type1.equalsIgnoreCase(type2))//this position doesn't match
                return false;
        }
        
        //if we get here they all match
        return true;
    }
    
    
    
    boolean canPrune(RCTuple cand, RCTuple comp, CheckSumType checkSumType){
        //see if competitive pruning is valid for the given candidate and competitor
        
        double checkSum = searchSpace.emat.getInternalEnergy(cand);
        
        checkSum -= searchSpace.emat.getInternalEnergy(comp);
        
        if(checkSumType == CheckSumType.GOLDSTEIN ){            
            
            for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
                if(!cand.pos.contains(pos)){
                    
                    double bestInteraction = Double.POSITIVE_INFINITY;
                    for(int rc: searchSpace.pruneMat.unprunedRCsAtPos(pos)){
                        
                        double diffBound = minInteraction(cand,pos,rc) - maxInteraction(comp,pos,rc);
                        
                        bestInteraction = Math.min(bestInteraction,
                               diffBound );
                    }
                    checkSum += bestInteraction;
                }
            }
        }
        else {
            throw new RuntimeException("ERROR: Not supporting indirect and conf-splitting pruning yet...");
        }
        
        return (checkSum>pruningInterval);
    }
    
    
    
    boolean canPrune(RCTuple cand, CheckSumType checkSumType){
        //try to prune cand non-competitively
                
        double checkSum = searchSpace.emat.getInternalEnergy(cand);
                
        if(checkSumType == CheckSumType.BOUNDS ){
            
            //this is like the lower bound in the A* ConfTree
            //break up the full energy into contributions associated with different res
            for(int level=0; level<searchSpace.confSpace.numPos; level++){
                double resContribLB = Double.POSITIVE_INFINITY;//lower bound on contribution of this residue
                //resContribLB will be the minimum_{rc} of the lower bound assuming rc assigned to this level
                if(!cand.pos.contains(level)){//level not fully defined
                    for ( int rc : searchSpace.pruneMat.unprunedRCsAtPos(level) ) {//cache this?
                        resContribLB = Math.min( resContribLB, RCContributionLB(level,rc,cand) );
                    }
                }
                
                checkSum += resContribLB;
            }
            
        }
        else {
            throw new RuntimeException("ERROR: Unrecognized checksum type for non-competitive pruning: "+checkSumType.name());
        }
        
        return ( checkSum > boundsThreshold+pruningInterval );
    }
    
    
    
    //For Bounds pruning.  Based on ConfTree.
    double RCContributionLB(int level, int rc, RCTuple definedTuple){
        //Provide a lower bound on what the given rc at the given level can contribute to the energy
        //assume definedTuple
        
        double rcContrib = 0;
        
        //for this kind of lower bound, we need to split up the energy into the defined-tuple energy
        //plus "contributions" for each undefined residue
        //so we'll say the "contribution" consists of any interactions that include that residue
        //but do not include higher-numbered undefined residues
        for(int level2=0; level2<level; level2++){
            
            if(definedTuple.pos.contains(level2) || level2<level){//lower-numbered or defined residues
                
                double levelBestE = Double.POSITIVE_INFINITY;//best pairwise energy
                
                ArrayList<Integer> allowedRCs = null;;
                if(definedTuple.pos.contains(level2)){
                    int index = definedTuple.pos.indexOf(level2);
                    int definedRC = definedTuple.RCs.get(index);
                    allowedRCs = new ArrayList<>();
                    allowedRCs.add(definedRC);
                }
                else
                    allowedRCs = searchSpace.pruneMat.unprunedRCsAtPos(level2);
                
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
    
    
    
    //For the rigid pruning with a pairwise matrix that is currently implemented,
    //min and maxInteraction are the same: just the sum of pairwise interactions
    //between (pos,rc) and all the RCs in RCTup
    double minInteraction(RCTuple RCTup, int pos, int rc){
        double E = 0;
        
        for(int indexInTup=0; indexInTup<RCTup.pos.size(); indexInTup++){
            int pos2 = RCTup.pos.get(indexInTup);
            int rc2 = RCTup.RCs.get(indexInTup);
            double pairwiseE = searchSpace.emat.getPairwise(pos, rc, pos2, rc2);
            E += pairwiseE;
        }
        
        return E;
    }
    
    
    double maxInteraction(RCTuple RCTup, int pos, int rc){
        return minInteraction(RCTup,pos,rc);
    }
    
    
    
    ArrayList<RCTuple> enumerateAllCandidates(int numPosInTuple){
        //enumerate all candidate tuples for pruning
        //basically all unpruned tuples of the specified size, at any position
        ArrayList<ArrayList<Integer>> posTupleCand = positionTuples(numPosInTuple);
        
        ArrayList<RCTuple> allCandidates = new ArrayList<>();
        
        for(ArrayList<Integer> posTuple : posTupleCand){
            ArrayList<RCTuple> candidatesAtPos = searchSpace.pruneMat.unprunedRCTuplesAtPos(posTuple);
            allCandidates.addAll(candidatesAtPos);
        }
        
        return allCandidates;
    }
    
    
    ArrayList<ArrayList<Integer>> positionTuples(int numPosInTuple){
        //get all possible sets of positions, of the specified size
        //(sets represented as tuples in ascending order, and must have all distinct positions)
        //
        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        int numPosTotal = searchSpace.confSpace.numPos;
        
        if(numPosInTuple==1){
            for(int pos=0; pos<numPosTotal; pos++){
                ArrayList<Integer> singleton = new ArrayList<>();
                singleton.add(pos);
                ans.add(singleton);
            }
        }
        else {
            ArrayList<ArrayList<Integer>> reducedTups = positionTuples(numPosInTuple-1);
            
            for(ArrayList<Integer> redTup : reducedTups){
                int lastPosInTup = redTup.get(numPosInTuple-2);
                
                for(int newPos=lastPosInTup+1; newPos<numPosTotal; newPos++){
                    ArrayList<Integer> fullTup = (ArrayList<Integer>)redTup.clone();
                    fullTup.add(newPos);
                    ans.add(fullTup);
                }
            }
        }
        
        return ans;
    }
    
    
    
    
}