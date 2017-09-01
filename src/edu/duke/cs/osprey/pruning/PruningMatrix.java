/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.*;

/**
 *
 * @author mhall44
 */
public class PruningMatrix extends TupleMatrixBoolean {
	
	private static final long serialVersionUID = 1212622649775905467L;
	
    //similar to energy matrix, but indicates what RCs and tuples of RCs are pruned
    //pruning indicated by true boolean
    //a conformation is pruned if it contains any pruned RC or tuple
    
    //private HigherTupleFinder[][][][] higherTerms;//look up higher terms by pair
    
    //maybe separate intra too?
    
	protected PruningMatrix(){//no allocation (for overriding by UpdatedPruningMatrix)
		super();
    }
	
	public PruningMatrix(PruningMatrix other) {
		super(other);
	}
        
    public PruningMatrix(ConfSpace cSpace, double pruningInterval){
        super(cSpace, pruningInterval, false);//higher tuples are unpruned unless otherwise indicated
        
        //We'll want to initialize everything to be unpruned, because this will be looked up during pruning
        //currently all entries in oneBody and pairwise are null
        fill(false);
    }

	public PruningMatrix(SimpleConfSpace confSpace) {
    	super(confSpace, 0, false);

    	// start with everything unpruned
		fill(false);
	}
    
    public PruningMatrix(int numPos, int[] numAllowedAtPos, double pruningInterval) {
    	super(numPos, numAllowedAtPos, pruningInterval, false);
    }
    
    public void unprunedRCsAtPos(ArrayList<Integer> out, int pos) {
    	out.clear();
    	int numRCs = getNumConfAtPos(pos);
		for (int index=0; index<numRCs; index++) {
			if(!getOneBody(pos,index))
				out.add(index);
		}
    }
    
    
    public ArrayList<Integer> unprunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> out = new ArrayList<>();
        unprunedRCsAtPos(out, pos);
        return out;
    }
    

    public void prunedRCsAtPos(ArrayList<Integer> out, int pos) {
    	out.clear();
    	int numRCs = getNumConfAtPos(pos);
		for (int index=0; index<numRCs; index++) {
			if(getOneBody(pos,index))
				out.add(index);
		}
    }
    
    
    public ArrayList<Integer> prunedRCsAtPos(int pos){
        //which RCs at the given position are unpruned?
        //Return index of the RCs within the position
        ArrayList<Integer> out = new ArrayList<>();
        prunedRCsAtPos(out, pos);
        return out;
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
    
    
    public boolean isPruned(RCTuple tup) {
    	
    	// OPTIMIZATION: this function gets hit a lot
    	// so even pedantic optimizations can have a noticeable impact
    	
    	// copy some references to stack
    	ArrayList<Integer> tuppos = tup.pos;
    	ArrayList<Integer> tupRCs = tup.RCs;
    	
    	// OPTIMIZATION: skipping even the check for higher order terms
    	// improves CPU cache performance a lot when we don't actually have any terms to use
    	boolean hasHigherOrderTerms = hasHigherOrderTerms();
    	
        //can be prune per se, or check if some singles in it are pruned, or pairs, etc.
    	
    	// OPTIMIZATION: checking singles separately from pairs improves CPU cache performance
    	// since we're memory-bound for this workload, trading the extra CPU instructions for
    	// fewer cache misses has a noticeable impact on performance
    	
    	// check singles
    	int numTupPos = tuppos.size();
        for (int i1=0; i1<numTupPos; i1++) {
            int pos1 = tuppos.get(i1);
            int rc1 = tupRCs.get(i1);
            
            if (getOneBody(pos1, rc1)) {
                return true;
            }
        }
            
        // check pairs
        for (int i1=0; i1<numTupPos; i1++) {
            int pos1 = tuppos.get(i1);
            int rc1 = tupRCs.get(i1);
            
            for (int i2=0; i2<i1; i2++) {
                int pos2 = tuppos.get(i2);
                int rc2 = tupRCs.get(i2);
            
                if (getPairwise(pos1, rc1, pos2, rc2)) {
                    return true;
                }
                
                if (hasHigherOrderTerms) {
					HigherTupleFinder<Boolean> htf = getHigherOrderTerms(pos1, rc1, pos2, rc2);
					if (htf != null) {
						if (isPrunedHigherOrder(tup, i2, htf)) {
							return true;
						}
					}
                }
            }
        }
        
        return false;
    }
    
    
    public boolean isPrunedHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Boolean> htf){
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

    public void pruneSingle(int pos1, int rc1) {
    	setOneBody(pos1, rc1, true);
    	prunePairsFromSingle(pos1, rc1);
	}

	public void prunePair(int pos1, int rc1, int pos2, int rc2) {
    	setPairwise(pos1, rc1, pos2, rc2, true);
	}

	public void prunePairsFromSingles() {
		int n = getNumPos();
		for (int pos1=0; pos1<n; pos1++) {
			int n1 = getNumConfAtPos(pos1);
			for (int rc1=0; rc1<n1; rc1++) {

				if (getOneBody(pos1, rc1)) {
					prunePairsFromSingle(pos1, rc1);
				}
			}
		}
	}

	public void prunePairsFromSingle(int pos1, int rc1) {
		int n = getNumPos();
		for (int pos2=0; pos2<n; pos2++) {

			if (pos1 == pos2) {
				continue;
			}
			int n2 = getNumConfAtPos(pos2);
			for (int rc2=0; rc2<n2; rc2++) {
				setPairwise(pos1, rc1, pos2, rc2, true);
			}
		}
	}

	public boolean isSinglePruned(int pos, int rc) {
    	return getOneBody(pos, rc);
	}

	public boolean isPairPruned(int pos1, int rc1, int pos2, int rc2) {
    	return getPairwise(pos1, rc1, pos2, rc2)
			|| isSinglePruned(pos1, rc1)
			|| isSinglePruned(pos2, rc2);
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
