/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author mhall44
 */
public class PruningMatrix extends TupleMatrix<Boolean> {
    //similar to energy matrix, but indicates what RCs and tuples of RCs are pruned
    //pruning indicated by true boolean
    //a conformation is pruned if it contains any pruned RC or tuple
    
    //private HigherTupleFinder[][][][] higherTerms;//look up higher terms by pair
    
    //maybe separate intra too?
    
        
    public PruningMatrix(ConfSpace cSpace, double pruningInterval){
        super(cSpace, pruningInterval);
        
        //We'll want to initialize everything to be unpruned, because this will be looked up during pruning
        //currently all entries in oneBody and pairwise are null
        for(ArrayList<Boolean> oneBodyAtPos : oneBody){
            Collections.fill(oneBodyAtPos, false);
        }
        
        for(ArrayList<ArrayList<ArrayList<Boolean>>> pairwiseAtPos : pairwise){
            for(ArrayList<ArrayList<Boolean>> pairwiseAt2Pos : pairwiseAtPos){
                for(ArrayList<Boolean> pairwiseAtRC : pairwiseAt2Pos){
                    Collections.fill(pairwiseAtRC, false);
                }
            }
        }
    }
    
    
    public ArrayList<Integer> unprunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> ans = new ArrayList<>();
        int numRCs = numRCsAtPos(pos);
         
        for(int index=0; index<numRCs; index++){
            if(!getOneBody(pos,index))
                ans.add(index);
        }
        
        return ans;
    }
    
    
    ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos){
        //get a list of unpruned RCTuples with the given positions
        //this method tests a few things more than once, so it could be sped up if needed, but it is convenient
        
        int numPos = pos.size();
        ArrayList<RCTuple> unpruned = new ArrayList<>();
        
        if(numPos==1){
            int posNum = pos.get(0);
            for(int rc=0; rc<oneBody.get(posNum).size(); rc++){
                if(!getOneBody(posNum,rc))
                    unpruned.add(new RCTuple(posNum,rc));
            }
        }
        else {
            //get unpruned tuples of RCs at all but the last position
            //then see what RCs at the last position we can add
            ArrayList<Integer> posReduced = (ArrayList<Integer>)pos.clone();
            posReduced.remove(numPos-1);
            
            ArrayList<RCTuple> tupsReduced = unprunedRCTuplesAtPos(posReduced);
            
            int lastPos = pos.get(numPos-1);
            
            for(int rc=0; rc<oneBody.get(lastPos).size(); rc++){
                if(!getOneBody(lastPos,rc)){
                    for(RCTuple reducedTup : tupsReduced){//try to combine into an unpruned RC
                        
                        ArrayList<Integer> fullRCList = (ArrayList<Integer>)reducedTup.RCs.clone();
                        fullRCList.add(rc);
                        
                        RCTuple fullTup = new RCTuple(pos,fullRCList);
                        if(!isPruned(fullTup))
                            unpruned.add(fullTup);
                    }
                }
            }
        }
        
        return unpruned;     
    }
    
    
    public boolean isPruned(RCTuple tup){
        //can be prune per se, or check if some singles in it are pruned, or pairs, etc.
        for(int indexInTup=0; indexInTup<tup.pos.size(); indexInTup++){
            int pos1 = tup.pos.get(indexInTup);
            int rc1 = tup.RCs.get(indexInTup);
            
            if(getOneBody(pos1,rc1))//rc1 pruned by itself
                return true;
            
            //check if any pairs pruned
            for(int index2=0; index2<indexInTup; index2++){
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
            
                if(getPairwise(pos1,rc1,pos2,rc2))
                    return true;
            }
            
            //check higher-order here...
        }
        return false;
    }
    
    
    public void markAsPruned(RCTuple tup){
        int tupSize = tup.pos.size();
        if(tupSize==1)
            setOneBody(tup.pos.get(0), tup.RCs.get(0), true);
        else if(tupSize==2)
            setPairwise(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), true);
        else
            throw new RuntimeException("ERROR: Can't record pruned tuples of size "+tupSize+" yet");
    }
    
    
    public int countPrunedRCs(){
        //how many RCs are pruned overall?
        int count = 0;
        for(ArrayList<Boolean> posPruned : oneBody){
            for(boolean b : posPruned){
                if(b)
                    count++;
            }
        }
        return count;
    }
    
    
    public int countPrunedPairs(){
        //how many pairs are pruned overall?
        int count = 0;
        for(ArrayList<ArrayList<ArrayList<Boolean>>> posPruned : pairwise){
            for(ArrayList<ArrayList<Boolean>> posPruned2 : posPruned){
                for(ArrayList<Boolean> RCPruned1 : posPruned2){
                    for(boolean b : RCPruned1){
                        if(b)
                            count++;
                    }
                }
            }
        }
        return count;
    }
    
    /*boolean isPruned(RC rc){
        //look up 1-body
        return getOneBody(pos,rcNum);
    }*/
    
    
    
    
}