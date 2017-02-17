/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import java.util.ArrayList;

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
        
        ArrayList<RCTuple> pairList = tupEnum.enumerateUnprunedTuples(2);
        if(pairList.isEmpty())//no pairs...indicates <2 flexible positions (or everything pruned, which will be detected)
            pairList.addAll(tupEnum.enumerateUnprunedTuples(1));//so add singles
        
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
        
        return expander.calcExpansion(tupList);
    }
            
}