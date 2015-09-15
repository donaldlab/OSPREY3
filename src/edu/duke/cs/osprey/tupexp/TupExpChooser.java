/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.Pruner;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * Chooses what tuple expansion to make for a particular problem (what tuples to include),
 * uses a TupleExpander to calculate it
 * 
 * Can use a pairwise expansion, or can look at how changes at certain positions affect 
 * other local energies, to get a comprehensive higher-order expansion
 * 
 * @author mhall44
 */
public class TupExpChooser {
    
    TupleExpander expander;
    TupleEnumerator tupEnum;

    
    
    public TupExpChooser(TupleExpander expander, TupleEnumerator tupEnum) {
        this.expander = expander;
        this.tupEnum = tupEnum;
    }
    
    
    public double calcPairwiseExpansion(){
        //calculate a full pairwise expansion
        
        ArrayList<RCTuple> pairList = tupEnum.enumerateUnprunedTuples(2);/*new ArrayList<>();
        
        
        //put in pairwise terms.  (1-body would be redundant)
        for(int pos=0; pos<expander.numPos; pos++){
            for(int a=0; a<expander.numAllowed[pos]; a++){
                
                RCTuple rc = new RCTuple(pos,a);
                        
                
                //We'll add all pairs that aren't pruned (as singles or pairs)
                //this ensures that all tuples have reasonable (and thus sample-able) conf spaces
                //because we'll only be adding higher tuples that have samples later on
                if(!expander.isPruned(rc)){

                    for(int pos2=0; pos2<pos; pos2++){

                        for(int a2=0; a2<expander.numAllowed[pos2]; a2++){
                            
                            RCTuple rc2 = new RCTuple(pos2,a2);
                            
                            if(!expander.isPruned(rc2)){
                                
                                RCTuple pair= new RCTuple(pos,a,pos2,a2);
                                
                                if(!expander.isPruned(pair)){
                                    pairList.add(pair);
                                }
                            }
                        }
                    }
                }
            }
        }*/
        
        //now calculate the expansion for this list of pairs
        return expander.calcExpansion(pairList);
    }
    
    
    public void calcExpansionCliqueTriples(double numTupFactor){
        //Augment the current pairwise expansion in expander using triples.  
        //We add at most numTupFactor * the current number of tuples
        //We look at triples with either two or three of their constituent pairs being unusually strong
        //in the pairwise expansion.  

        int numTriplesPerType = (int) (numTupFactor * expander.tuples.size() / 2);//max num triples to consider of each type
        
        ArrayList<RCTuple> tupList = (ArrayList<RCTuple>) expander.tuples.clone();
        ArrayList<RCTuple> topTriples = tupEnum.clique2PairTriples(numTriplesPerType);
        tupList.addAll(topTriples);
        
        expander.calcExpansion(tupList);
    }
    
    
    
    public double calcExpansionResTriples(int numPartners){
        //we'll consider all triples at residue triples of special interest
        //(i.e., containing at least two pairs that are top interactions for at least one of the residues in the pair)
        
        
        if(tupEnum.getEmat() == null)//no base energy matrix...use pairwise expansion here
            tupEnum.setEmat( expander.getEnergyMatrix() );
        
        
        ArrayList<RCTuple> tupList = (ArrayList<RCTuple>) expander.tuples.clone();
        ArrayList<ArrayList<Integer>> topPositionTriples = tupEnum.topPositionTriples(numPartners);
        ArrayList<RCTuple> topTriples = tupEnum.enumerateUnprunedTuples(topPositionTriples);
        tupList.addAll(topTriples);
        
        
        
        //Now, use all the RC triples that involve at least two "strong" interactions
        //DEBUG!!
        /*Pruner pruner = new Pruner( ((ConfETupleExpander)expander).sp, false, 0, 0 );
        ArrayList<RCTuple> possibleTriples = pruner.enumerateAllCandidates(3);
                
        ArrayList<RCTuple> tupList = (ArrayList<RCTuple>) expander.tuples.clone();
        
        //go through possible triples
        for(RCTuple triple : possibleTriples) {
            int strongInteractionCount = 0;
            for(int i=0; i<3; i++){
                int p1 = triple.pos.get(i);
                int p2 = triple.pos.get((i+1)%3);
                if(strongInteraction[p1][p2])
                    strongInteractionCount++;
           }
           
           if(strongInteractionCount>=2)
               tupList.add(triple);
        }*/
        
        return expander.calcExpansion(tupList);
    }
            
    
    
    
    /*
    public void calcExpansion(double errorThresh){
        //approximate within the specified error, using interaction cutoffs
        
        double interactionCutoff = initEstCutoff(errorThresh,expander.numPos);
        //hmm if resid thresh is errorThresh, overall conf errors can be about sqrt(errorThresh),
        //which is split over maybe 4*numPos or so major interactions?  So sqrt(errorThresh)/(4*numPos)?
        //but try to high-ball it since we can reduce interaction cutoff if needed?  sqrt(errorThresh)/(2*numPos)?
           
        double error = Double.POSITIVE_INFINITY;
        
        System.out.println("Choosing tuple expansion by interaction cutoffs, based on error"
                + " threshold of "+errorThresh);
        
        
        while(error>errorThresh){
            
            ArrayList<RCTuple> tupleList = listInteractingTuples(interactionCutoff);
            
            System.out.println("Interaction cutoff: "+interactionCutoff);
            
            error = expander.calcExpansion(tupleList);//OR UPDATE?
            
            interactionCutoff /= 2;//we'll reduce until we get a good error
        }
        
        System.out.println("Tuple expansion reached satisfactory residual: "+error);
        
        //when we get here, the expansion is good enough
        //we'll always get here, because if the cutoff is small enough all interactions will be included
    }
    
    
    
    double initEstCutoff(double errorThresh, int numPos){
        return 0.5*Math.sqrt(errorThresh)/numPos;
    }
            
    
    //FOR HIGHER-ORDER WILL NEED TO MODIFY CONFTREE, PRUNER, AND E/PRUNE MATS
    //IN ADDITION TO TUP EXP STUF
    
    
    //These functions select tuples based on estimates of local interactions
    //The choices are validate with an actual tuple expansion, so just roughly getting it right is sufficient
    //we do need to be able to eventually include any significant interactions by lowering the interaction energy cutoff though
    
    //PROBLEM WITH APPROACH BELOW IS THAT MOST MODULATIONS OF TUPLE 1 BASE ENERGY BY TUPLE 2
    //CAN BE EXPLAINED IN TERMS OF TUPLES SMALLER THAN UNION(TUPLE 1, TUPLE 2)
    
    ArrayList<RCTuple> listInteractingTuples(double cutoff){
        //What tuples have interactions above some cutoff?
        //The "interactions" in a tuple (a,b,...,i,j,...) are differences in the local energy for i,j,...
        //due to changing RCs at a,b,... that are NOT accounted for by changing sub-combinations of a,b,...
        
        ArrayList<RCTuple> tupleList = new ArrayList<>();
        RCTuple nothing = new RCTuple();
        
        for(RCTuple baseTup : baseInteractingTuples()){
            
            if(expander.isPruned(baseTup))
                continue;
            
            HashSet<RCTuple> interactionsForBase = new HashSet<>();//list of RC tuples significantly involving the base tuple...
            //WATCH EQUALS
            new TupInteractionsEstimator().checkInteractions(baseTup);
            //will generate samples with tuple
            //if can change a,b,... in sample, then will do that
            //if not, will relax clashes involving a,b,... (preferably not also i,j,but yes if needed) and try again
            //(relaxed-clash min. can have its own function--not within purview of TupleExpander, involves E fcn
            //construction)
            
            if(enough){
                check AA split, if so within AA;
                recurse;
            }
        }
        
        return tupleList;
    }
    
    
    
    void addInteractingTuples(RCTuple curTup, ArrayList<RCTuple> tupleList, double cutoff){
        //figure out if curTup has enough interactions (based on cutoff) to include in tupleList
        
        
        //if we get here there are either no or sub-cutoff interactions in the base energy model (i.e., rigid energies real small)...

        checkTupleInteractions(curTup);
        
        //try AA-based splitting.  Possibly move on to all RCs
        
    }
    
    boolean checkTupleInteractions(RCTuple tup, double cutoff, boolean withinAA){
        
        //does tup have any significant interactions beyond what sub-tuples can offer?
        //if there are clashes, just drop those interactions consistently
        //MIGHT WANT A SPECIAL CLASS...THIS IS SORT OF SPECIFIC TO WHAT'S BEING EXPANDED

        if(curTup.hasBaseInteractions()) {
            if bigEnough) return yeah;
        }


        for(10 samples):
            assign tup
            TESampleSet.finish()
            foreach sub-tuple st that has local interactions:
                randomly change RCs for tup-st (withinAA if indicated)
    }
    */
}
