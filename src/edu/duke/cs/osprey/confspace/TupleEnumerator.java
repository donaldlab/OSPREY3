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

package edu.duke.cs.osprey.confspace;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * @author mhall44
 */
public class TupleEnumerator {
    //Makes lists of residue or RC tuples for pruning, tuple expansion, etc.
    
    PruningMatrix pruneMat;//for figuring out which are unpruned
    EnergyMatrix emat;//for figuring out which tuples interact strongly
    int numPosTotal;//total number of positions

    
    public TupleEnumerator(PruningMatrix pruneMat, EnergyMatrix emat, int numPosTotal) {
        this.pruneMat = pruneMat;
        this.emat = emat;
        this.numPosTotal = numPosTotal;
    }
    
    
    
    public ArrayList<RCTuple> enumerateUnprunedTuples(int numPosInTuple){
        //enumerate all unpruned tuples of the specified size, at any position
        //e.g. if numPosInTuple==2 then enumerate all unpruned pairs
        ArrayList<ArrayList<Integer>> posTupleCand = allPositionTuples(numPosInTuple);
        return enumerateUnprunedTuples( posTupleCand );
    }
    
    
    public ArrayList<RCTuple> enumerateUnprunedTuples( ArrayList<ArrayList<Integer>> posTupleCand  ){
        //Enumerate all unpruned RC tuples at the specified position tuples
        
        ArrayList<RCTuple> allCandidates = new ArrayList<>();
        
        for(ArrayList<Integer> posTuple : posTupleCand){
            ArrayList<RCTuple> candidatesAtPos = pruneMat.unprunedRCTuplesAtPos(posTuple);
            allCandidates.addAll(candidatesAtPos);
        }
        
        return allCandidates;
    }
    
    
    public ArrayList<ArrayList<Integer>> allPositionTuples(int numPosInTuple){
        //get all possible sets of positions, of the specified size
        //(sets represented as tuples in ascending order, and must have all distinct positions)

        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        
        if(numPosInTuple==1){
            for(int pos=0; pos<numPosTotal; pos++){
                ArrayList<Integer> singleton = new ArrayList<>();
                singleton.add(pos);
                ans.add(singleton);
            }
        }
        else {
            ArrayList<ArrayList<Integer>> reducedTups = allPositionTuples(numPosInTuple-1);
            
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
    
    
    public ArrayList<ArrayList<Integer>> topPositionTriples(int numPartners){
        //Get only the strongest-interacting position triples
        
        //enumerate all interactions first...  (pretty quick because not looping over RCs)
        ArrayList<ArrayList<Integer>> allPositionTriples = allPositionTuples(3);
        ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
        
        boolean strongInteraction[][] = strongInteractingPairs(numPartners);
        
        //pick out the strong ones
        for(ArrayList<Integer> posTriple : allPositionTriples) {
            int strongInteractionCount = 0;
            for(int i=0; i<3; i++){
                int p1 = posTriple.get(i);
                int p2 = posTriple.get((i+1)%3);
                if(strongInteraction[p1][p2])
                    strongInteractionCount++;
           }
           
           if(strongInteractionCount>=2)
               ans.add(posTriple);
        }
        
        return ans;
    }
    
    
    public boolean[][] strongInteractingPairs(int numPartners){
        //For each position pos, the numPartners other positions that interact most strongly with pos
        //are counted as "strong" interactions
        //return a matrix of which pairwise interactions are "strong"
        //Using lower-bound energy matrix should capture these trends without possible noise from initial bad tup-exp fit...
        
        int numPos = emat.getNumPos();
        double strongestPairE[][] = emat.topPairwiseInteractions();
        
        //Next, use these to figure out what are the numPartners top interaction partners for each residue
        //make a matrix of "strong" interactions based on this
        boolean strongInteraction[][] = new boolean[numPos][numPos];
        
        for(int pos=0; pos<numPos; pos++){
            
            PriorityQueue<TupE> posTop = new PriorityQueue<>();
            
            for(int pos2=0; pos2<numPos; pos2++){
                if(pos!=pos2){
                    posTop.add( new TupE(new RCTuple(pos2,-1), strongestPairE[pos][pos2]) );
                    if(posTop.size()>numPartners)//throw out subpar interactions
                        posTop.poll();
                }
            }
            
            for(TupE top : posTop){
                int pos2 = top.tup.pos.get(0);
                strongInteraction[pos][pos2] = true;
                strongInteraction[pos2][pos] = true;
            }
        }
        
        return strongInteraction;
    }
    
    
    public ArrayList<RCTuple> clique2PairTriples(int numTriplesPerType){
        //enumerate RC triples based on containing either two or three strong pairwise interactions in emat
        //we take the best triples for each measure (numTriplesPerType of each)
        
        ArrayList<RCTuple> tupList = new ArrayList<>();
        
        PriorityQueue<TupE> topCliques = new PriorityQueue<>();//top triples in terms of all pairs strong
        PriorityQueue<TupE> top2Pair = new PriorityQueue<>();//top triples in terms of having 2 strong pairs
        
        ArrayList<RCTuple> possibleTriples = enumerateUnprunedTuples(3);//consider all unpruned triples
        
        //go through possible triples
        for(RCTuple triple : possibleTriples) {
            double[] pairAbsE = new double[3];
            double minAbsE = Double.POSITIVE_INFINITY;
            double max2PairE = Double.NEGATIVE_INFINITY;
            
            for(int i=0; i<3; i++){
                RCTuple pair = triple.subtractMember(i);
                pairAbsE[i] = Math.abs( emat.getPairwise(pair.pos.get(0), pair.RCs.get(0), 
                        pair.pos.get(1), pair.RCs.get(1)) );
                
                minAbsE = Math.min(minAbsE,pairAbsE[i]);
            }
            
            for(int i=0; i<3; i++)
                max2PairE = Math.max( max2PairE, Math.min(pairAbsE[(i+1)%3],pairAbsE[(i+2)%3]) );
            
            
            topCliques.add( new TupE(triple,minAbsE) );
            top2Pair.add( new TupE(triple,max2PairE) );
            
            //now throw out the weakest-interacting triples if we have too many
            if(topCliques.size()>numTriplesPerType)
                topCliques.poll();//CHECK DIRECTIONS
            if(top2Pair.size()>numTriplesPerType)
                top2Pair.poll();
        }
        
        
        for(TupE tupe : topCliques)
            tupList.add(tupe.tup);
        for(TupE tupe : top2Pair)
            tupList.add(tupe.tup);
        
        return tupList;
    }

    public EnergyMatrix getEmat() {
        return emat;
    }

    public void setEmat(EnergyMatrix emat) {
        this.emat = emat;
    }
    
    
    
    
}
