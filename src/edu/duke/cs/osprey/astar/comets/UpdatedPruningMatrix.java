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

package edu.duke.cs.osprey.astar.comets;

import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * Within a given state, each node in the COMETS tree can have a different set of unpruned 
 * RCs and pairs, but these will usually be small updates to the previous node
 * 
 * This "updated" matrix will behave as a pruning matrix for purposes of storing/checking
 * pruned RCs and pairs at a node, but without full storage
 * 
 * @author mhall44
 */
public class UpdatedPruningMatrix extends PruningMatrix {
    
    public PruningMatrix parent;//This matrix will have everything pruned in parent, plus updates
    
    ArrayList<TreeSet<Integer>> prunedRCUpdates = new ArrayList<>();
    //list of pruned RCs (just for this update) at each position

    
    ArrayList<ArrayList<TreeMap<Integer,TreeSet<Integer>>>> prunedPairUpdates = new ArrayList<>();
    //list of pruned RCs (just for this update) at each position
    
    
    public UpdatedPruningMatrix(PruningMatrix parent) {
        this.parent = parent;
        
        int numPos = parent.getNumPos();
        
        for(int pos=0; pos<numPos; pos++){
            prunedRCUpdates.add(new TreeSet<Integer>());
            prunedPairUpdates.add(new ArrayList<TreeMap<Integer,TreeSet<Integer>>>());
            
            for(int pos2=0; pos2<pos; pos2++){
                prunedPairUpdates.get(pos).add(new TreeMap<Integer,TreeSet<Integer>>());
            }
        }
    }
    
    
    
    @Override
    public void markAsPruned(RCTuple tup){
        //Store as update
        int tupNumPos = tup.pos.size();
        
        if(tupNumPos==1){
            int pos = tup.pos.get(0);
            int rc =  tup.RCs.get(0);
            prunedRCUpdates.get(pos).add(rc);
        }
        else if(tupNumPos==2){
            int pos1 = tup.pos.get(0);
            int pos2 = tup.pos.get(1);
            int rc1 =  tup.RCs.get(0);
            int rc2 =  tup.RCs.get(1);
            
            if(pos1<pos2){//need to store the pair in descending order of position
                pos2 = tup.pos.get(0);
                pos1 = tup.pos.get(1);
                rc2 =  tup.RCs.get(0);
                rc1 =  tup.RCs.get(1);
            }
            
            TreeMap<Integer,TreeSet<Integer>> pairs = prunedPairUpdates.get(pos1).get(pos2);
            
            if(!pairs.containsKey(rc1))//allocate the treeset for pairs involving rc1
                pairs.put(rc1, new TreeSet<Integer>());
            
            pairs.get(rc1).add(rc2);
        }
        else{
            throw new RuntimeException("ERROR: UpdatedPruningMatrix just stores updated"
                    + " singles and pairs pruning, can't store pruned tuple: "+tup.stringListing());
        }
    }
    
        
    @Override
    public Boolean getPairwise(int res1, int index1, int res2, int index2){
        //working with residue-specific RC indices directly.  
        
        if(parent.getPairwise(res1, index1, res2, index2))//first check parent
            return true;
        
        //also check updates
        if(res1>res2)
            return checkIntPair( prunedPairUpdates.get(res1).get(res2), index1, index2 );
        else
            return checkIntPair( prunedPairUpdates.get(res2).get(res1), index2, index1 );
    }
    
    
    private static boolean checkIntPair(TreeMap<Integer,TreeSet<Integer>> pairs, int rc1, int rc2){
        //Check if (rc1, rc2) is in the map (ordered tuple)
        if(pairs.containsKey(rc1)){
            if(pairs.get(rc1).contains(rc2))
                return true;
        }
        return false;
    }
    
    
    @Override
    public Boolean getOneBody(int res, int index){
        
        if(parent.getOneBody(res,index))//first check parent
            return true;
        
        //also check updates
        return prunedRCUpdates.get(res).contains(index);
    }
    
    
    //No updates for higher order
    @Override
    public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2){
        return parent.getHigherOrderTerms(res1, index1, res2, index2);
    }
        
    
    
    public int countUpdates(){
        //How many update RCs and pairs are there, put together?
        int count = 0;
        
        for(TreeSet<Integer> posUpdates : prunedRCUpdates){
            count += posUpdates.size();
        }
        
        for(ArrayList<TreeMap<Integer,TreeSet<Integer>>> posUpdates : prunedPairUpdates){
            for(TreeMap<Integer,TreeSet<Integer>> ppUpdates : posUpdates){
                for(TreeSet<Integer> pUpdates : ppUpdates.values()){
                    count += pUpdates.size();
                }
            }
        }
        
        return count;
    }
    
    
    @Override
    public int getNumConfAtPos(int pos){
        return parent.getNumConfAtPos(pos);
    }
    
    
    @Override
    public int getNumPos(){
        return parent.getNumPos();
    }
}
