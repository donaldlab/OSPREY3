/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;

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
    
    public PruningMatrix(){//no allocation (for overriding by UpdatedPruningMatrix)
        super(false);
    }
        
    public PruningMatrix(ConfSpace cSpace, double pruningInterval){
        super(cSpace, pruningInterval, false);//higher tuples are unpruned unless otherwise indicated
        
        //We'll want to initialize everything to be unpruned, because this will be looked up during pruning
        //currently all entries in oneBody and pairwise are null
        fill(false);
    }
    
    
    public ArrayList<Integer> unprunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> ans = new ArrayList<>();
        int numRCs = getNumConfAtPos(pos);
         
        for(int index=0; index<numRCs; index++){
            if(!getOneBody(pos,index))
                ans.add(index);
        }
        
        return ans;
    }
    
    
    public ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos){
        //get a list of unpruned RCTuples with the given positions
        //this method tests a few things more than once, so it could be sped up if needed, but it is convenient
        
        int numPos = pos.size();
        ArrayList<RCTuple> unpruned = new ArrayList<>();
        
        if(numPos==1){
            int posNum = pos.get(0);
            for(int rc=0; rc<getNumConfAtPos(posNum); rc++){
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
            
            for(int rc=0; rc<getNumConfAtPos(lastPos); rc++){
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
                
                HigherTupleFinder<Boolean> htf = getHigherOrderTerms(pos1,rc1,pos2,rc2);
                if(htf != null){
                    if(isPrunedHigherOrder(tup,index2,htf))
                        return true;
                }
            }
            
        }
        return false;
    }
    
    
    boolean isPrunedHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Boolean> htf){
        //Checks if tup is pruned based on interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        
        for(int ipos : interactingPos){
            
            //see if ipos is in tup with index < curIndex
            int iposIndex = -1;
            for(int ind=0; ind<curIndex; ind++){
                if(tup.pos.get(ind)==ipos){
                    iposIndex = ind;
                    break;
                }
            }

            if(iposIndex > -1){//ipos interactions need to be counted
                int iposRC = tup.RCs.get(iposIndex);
                if( htf.getInteraction(ipos, iposRC) )//sub-tuple plus (ipos,iposRC) is pruned
                    return true;
                
                //see if need to go up to highers order again...
                HigherTupleFinder htf2 = htf.getHigherInteractions(ipos,iposRC);
                if(htf2!=null){
                    if( isPrunedHigherOrder(tup,iposIndex,htf2) )
                        return true;
                }
            }
        }
        
        //if we get here, not pruned
        return false;
    }
    
    
    
    public void markAsPruned(RCTuple tup){
        setTupleValue(tup, true);
        /*
        int tupSize = tup.pos.size();
        if(tupSize==1)
            setOneBody(tup.pos.get(0), tup.RCs.get(0), true);
        else if(tupSize==2)
            setPairwise(tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), true);
        else
        */
    }
    
    
    public int countPrunedRCs(){
        //how many RCs are pruned overall?
        int count = 0;
        int numPos = getNumPos();
        for (int res1=0; res1<numPos; res1++) {
        	int m1 = getNumConfAtPos(res1);
        	for (int i1=0; i1<m1; i1++) {
        		if (getOneBody(res1, i1) == true) {
        			count++;
        		}
        	}
        }
        return count;
    }
    
    
    public int countPrunedPairs(){
        //how many pairs are pruned overall?
        int count = 0;
        int numPos = getNumPos();
        for (int res1=0; res1<numPos; res1++) {
        	int m1 = getNumConfAtPos(res1);
        	for (int i1=0; i1<m1; i1++) {
        		for (int res2=0; res2<res1; res2++) {
        			int m2 = getNumConfAtPos(res2);
        			for (int i2=0; i2<m2; i2++) {
        				if (getPairwise(res1, i1, res2, i2) == true) {
        					count++;
        				}
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