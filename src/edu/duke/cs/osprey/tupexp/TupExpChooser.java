/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
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
