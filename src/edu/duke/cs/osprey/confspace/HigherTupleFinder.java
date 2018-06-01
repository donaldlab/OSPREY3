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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * This is an object associated to a particular tuple t (pair or higher) in a TupleMatrix
 * It lists any higher-tuples that include t and that have non-trivial values
 * ("trivial" meaning 0 for energies, unpruned for pruning)
 * Tuples of order size(t)+1 are held explicitly; we recursively list HigherTupleFinder
 * objects for higher orders if there are even bigger tuples
 * 
 * @author mhall44
 */
public class HigherTupleFinder<T> implements Serializable {
    
    private ArrayList<Integer> interactingPos = new ArrayList<>();
    //flexible positions interacting with t in higher tuples
    
    //For sparseness purposes, the following lists are indexed by interacting position,
    //not by the usual position number
    //They map RC numbers (at the specified position) to the value of interest
    private ArrayList<TreeMap<Integer,T>> interactions = new ArrayList<>();
    private ArrayList<TreeMap<Integer,HigherTupleFinder<T>>> higher = new ArrayList<>();//next order up...
    
    //Note: the RC number could potentially be a code (e.g., -2) indicating a set of RCs,
    //not just one RC.
    //if so the TupleMatrix should store some indication of what these codes are
    
    T defaultInteraction;//If no interaction included explicitly, this is the default value
    //(e.g., 0 for energies, false for pruning)
    
    public HigherTupleFinder(T defaultInteraction){
        //generally start empty
        this.defaultInteraction = defaultInteraction;
    }

    public ArrayList<Integer> getInteractingPos() {
        return interactingPos;
    }
    
    public T getInteraction(int posNum, int RCNum) {
        //get interaction of this tuple with the RC (posNum,RCNum)
        
        for(int interactionNum=0; interactionNum<interactingPos.size(); interactionNum++){
            //there should be only a few interacting positions so we can loop over them quickly
            
            if(interactingPos.get(interactionNum) == posNum){
                
                if(interactions.get(interactionNum).containsKey(RCNum))
                    return interactions.get(interactionNum).get(RCNum);
                else
                    return defaultInteraction;
            }
        }
        
        //if we get here no interaction (i.e., just default)
        return defaultInteraction;
    }
    
    public HigherTupleFinder<T> getHigherInteractions(int posNum, int RCNum) {
        //get higher-order interactions involving the super-tuple 
        //consisting of this tuple plus the RC (posNum,RCNum)
        
        for(int interactionNum=0; interactionNum<interactingPos.size(); interactionNum++){
            //there should be only a few interacting positions so we can loop over them quickly
            
            if(interactingPos.get(interactionNum) == posNum){
                
                if(higher.get(interactionNum).containsKey(RCNum))
                    return higher.get(interactionNum).get(RCNum);
                else
                    return null;
            }
        }
        
        //if we get here no higher interaction
        return null;
    }
    
    
    
    public void setInteraction(RCTuple tup, T val){
        //set the interaction of this tuple with tup to the given value
        
        if(tup.pos.size()==1){//store interaction directly in this HigherTupleFinder
            int pos = tup.pos.get(0);
            int rc = tup.RCs.get(0);
            
            int posIndex = getPosIndex(pos);
            
            interactions.get(posIndex).put(rc, val);
        }
        else {//kick it up to the next level.  Make sure all sub-tuples of tup know about this interaction.  
            
            for(int index=0; index<tup.pos.size(); index++){
                int pos = tup.pos.get(index);
                int rc = tup.RCs.get(index);
                
                int posIndex = getPosIndex(pos);
                
                RCTuple subTup = tup.subtractMember(index);
                
                HigherTupleFinder<T> nextHTF = higher.get(posIndex).get(rc);
                
                if(nextHTF == null){//allocate next-level HTF if not currently existent
                    nextHTF = new HigherTupleFinder<>(defaultInteraction);
                    higher.get(posIndex).put(rc, nextHTF);
                }
                
                nextHTF.setInteraction(subTup, val);
            }
        }
    }
    
    
    private int getPosIndex(int pos){
        //find the index in interactingPos of pos
        //Create one if non-existent
        //This function is called when we need to set new interactions
        
        for(int interactionNum=0; interactionNum<interactingPos.size(); interactionNum++){
            //there should be only a few interacting positions so we can loop over them quickly
            
            if(interactingPos.get(interactionNum) == pos)
                return interactionNum;
        }
        
        //not found...create
        interactingPos.add(pos);
        interactions.add(new TreeMap<Integer,T>());
        higher.add(new TreeMap<Integer,HigherTupleFinder<T>>());
        return interactingPos.size()-1;
    }
    
    
    public ArrayList<RCTuple> listInteractionsWithValue(T val){
        //list tuples with the given interaction value
        //stored here
        
        ArrayList<RCTuple> ans = new ArrayList<>();
        
        //search recursively, recording tuples in decsending order of position
        recordInteractionsWithValue(val,ans,Integer.MAX_VALUE);
        
        return ans;
    }
    
    
    private void recordInteractionsWithValue(T val, ArrayList<RCTuple> tupList, int maxPos){
        //Find tuples with all pos numbers less than maxPos that have interactions value val
        //add them to tupList
        
        for(int interactionNum=0; interactionNum<interactingPos.size(); interactionNum++){
            
            int pos = interactingPos.get(interactionNum);
            
            if(pos<maxPos){
                
                for(int rc : interactions.get(interactionNum).keySet()){
                    
                    if(interactions.get(interactionNum).get(rc) == val){
                        tupList.add( new RCTuple(pos,rc) );
                    }
                }
                
                for(int rc : higher.get(interactionNum).keySet()){
                    
                    ArrayList<RCTuple> subTupList = new ArrayList<>();
                    higher.get(interactionNum).get(rc).recordInteractionsWithValue(val, subTupList, pos);
                    //the HigherTupleFinders in higher don't know about (pos,rc), so add that
                    
                    for(RCTuple subTup : subTupList){
                        subTup.pos.add(pos);
                        subTup.RCs.add(rc);
                        tupList.add(subTup);
                    }
                }
            }
        }
    }
    
    
}
