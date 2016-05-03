/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */
public class TupleMatrix<T> implements Serializable {
	
	private static final long serialVersionUID = 845854459137269739L;
	
    //We will need "matrices" of quantities defined
    //for example, the energy matrix (T=Double) stores single, pairwise, and higher-order energies
    //and we'll also have pruning (T=Boolean) and EPIC (T=EPoly) matrices
    //we'll store things as ArrayLists to make it easier to merge residues, partition RCs, etc.
    //and also to facilitate generics
	
	private int numPos; // eg residues
	private int[] numConfAtPos; // eg RCs at each residue
    
    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
	
	private ArrayList<T> oneBody; // indices: res1, RC1
	private ArrayList<T> pairwise; // indices: res1, res2, RC1, RC2 where res1>res2
    private ArrayList<HigherTupleFinder<T>> higherTerms; // indices: same as pairwise, can be null if no interactions
    
    // index the arrays
    private int[] oneBodyOffsets;
    private int[] pairwiseOffsets;
    
    //maybe separate intra too?
    
    
    //The above all use RCs that are indexed by residues
    //the following lists indicate exactly what RCs those are
    //private because changing it may require re-indexing everything else
    //METHODS TO ADD RCs?
    
    //private ArrayList<ArrayList<RC>> RCList;
    
    /*
    //reverse lookup: for each residue, a lookup first by AA type, then by 
    //Return the residue-based RC number
    ArrayList<TreeMap<String,TreeMap<Integer,
    */
    
    
    private double pruningInterval;//This matrix needs to hold entries for all RCs
    //that cannot be pruned with the specified pruning interval (Ew + Ival)
    //i.e. the matrix must describe all conformations within pruningInterval 
    //of the lowest pairwise lower bound
    
    private T defaultHigherInteraction;//We only mark sparse higher interactions;
    //if unmarked we assume this value (e.g., 0 for energy, false for pruning)
    
    
    protected TupleMatrix(T defaultHigherInteraction){
        //no allocation (for overriding)
        this.defaultHigherInteraction = defaultHigherInteraction;
    }
    
    
    public TupleMatrix(ConfSpace cSpace, double pruningInterval, T defaultHigherInteraction){
        //allocate the matrix based on the provided conformational space
        init(cSpace.numPos, cSpace.getNumRCsAtPos(), pruningInterval, defaultHigherInteraction);
    }
    
    
    public TupleMatrix(int numPos, int[] numAllowedAtPos, double pruningInterval, T defaultHigherInteraction){
        //allocate the matrix based on the provided conformational space size
        //also specify what pruningInterval it's valid up to
        init(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }
    
    
    private void init(int numPos, int[] numConfAtPos, double pruningInterval, T defaultHigherInteraction) {
        
    	this.numPos = numPos;
    	this.numConfAtPos = numConfAtPos;
        this.pruningInterval = pruningInterval;
        this.defaultHigherInteraction = defaultHigherInteraction;
        
        // compute the indices and allocate space
        int offset = 0;
        
        // first one-body offsets
        oneBodyOffsets = new int[numPos];
        for (int res1=0; res1<numPos; res1++) {
        	oneBodyOffsets[res1] = offset;
        	offset += numConfAtPos[res1];
        }
        
        // allocate space
        oneBody = new ArrayList<>(offset);
        for (int i=0; i<offset; i++) {
        	oneBody.add(null);
        }
        
        // then pairwise offsets
        pairwiseOffsets = new int[numPos*(numPos - 1)/2];
        offset = 0;
        int pairwiseIndex = 0;
        for (int res1=0; res1<numPos; res1++) {
        	for (int res2=0; res2<res1; res2++) {
        		pairwiseOffsets[pairwiseIndex++] = offset;
        		offset += numConfAtPos[res1]*numConfAtPos[res2];
        	}
        }
        assert (pairwiseIndex == pairwiseOffsets.length);
        
        // allocate space
        pairwise = new ArrayList<>(offset);
        // TODO: use lazy allocation for higher terms
        higherTerms = new ArrayList<>(offset);
        for (int i=0; i<offset; i++) {
        	pairwise.add(null);
        	higherTerms.add(null);
        }
    }
    
    public void fill(T val) {
		for (int res1=0; res1<numPos; res1++) {
			int m1 = numConfAtPos[res1];
			for (int i1=0; i1<m1; i1++) {
				setOneBody(res1, i1, val);
				for (int res2=0; res2<res1; res2++) {
					int m2 = numConfAtPos[res2];
					for (int i2=0; i2<m2; i2++) {
						setPairwise(res1, i1, res2, i2, val);
					}
				}
			}
		}
    }
    
    public void fill(Iterator<T> val) {
		for (int res1=0; res1<numPos; res1++) {
			int m1 = numConfAtPos[res1];
			for (int i1=0; i1<m1; i1++) {
				setOneBody(res1, i1, val.next());
				for (int res2=0; res2<res1; res2++) {
					int m2 = numConfAtPos[res2];
					for (int i2=0; i2<m2; i2++) {
						setPairwise(res1, i1, res2, i2, val.next());
					}
				}
			}
		}
    }
    
    public double getPruningInterval() {
        return pruningInterval;
    }
    
    public int getNumPos() {
    	return numPos;
    }
    
    public int getNumConfAtPos(int pos) {
    	return numConfAtPos[pos];
    }
    
    private int getOneBodyIndex(int res, int conf) {
    	return oneBodyOffsets[res] + conf;
    }
    
    private int getPairwiseIndexNoCheck(int res1, int res2) {
    	return res1*(res1 - 1)/2 + res2;
    }
    
    private int getPairwiseIndex(int res1, int res2) {
    	
    	// res2 should be strictly less than res1
    	if (res2 > res1) {
    		int swap = res1;
    		res1 = res2;
    		res2 = swap;
    	} else if (res1 == res2) {
    		throw new Error("Can't pair residue " + res1 + " with itself");
    	}
    	
    	return getPairwiseIndexNoCheck(res1, res2);
    }
    
    private int getPairwiseIndex(int res1, int conf1, int res2, int conf2) {
    	
    	// res2 should be strictly less than res1
    	if (res2 > res1) {
    		int swap = res1;
    		res1 = res2;
    		res2 = swap;
    		swap = conf1;
    		conf1 = conf2;
    		conf2 = swap;
    	} else if (res1 == res2) {
    		throw new Error("Can't pair residue " + res1 + " with itself");
    	}
    	
    	return pairwiseOffsets[getPairwiseIndexNoCheck(res1, res2)] + numConfAtPos[res2]*conf1 + conf2;
    }
    
    public T getOneBody(int res, int conf) {
    	return oneBody.get(getOneBodyIndex(res, conf));
    }
    
    public void setOneBody(int res, int conf, T val) {
    	oneBody.set(getOneBodyIndex(res, conf), val);
    }
    
    public void setOneBody(int res, ArrayList<T> val) {
    	int n = numConfAtPos[res];
    	int firstIndex = oneBodyOffsets[res];
    	for (int i=0; i<n; i++) {
    		oneBody.set(firstIndex + i, val.get(i));
    	}
    }
    
    public T getPairwise(int res1, int conf1, int res2, int conf2) {
    	return pairwise.get(getPairwiseIndex(res1, conf1, res2, conf2));
    }
    
    public void setPairwise(int res1, int conf1, int res2, int conf2, T val) {
    	pairwise.set(getPairwiseIndex(res1, conf1, res2, conf2), val);
    }
    
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<T>> val) {
    	
    	// res2 should be strictly less than res1
    	if (res2 > res1) {
    		int swap = res1;
    		res1 = res2;
    		res2 = swap;
    	} else if (res1 == res2) {
    		throw new Error("Can't pair residue " + res1 + " with itself");
    	}
    	
    	int n1 = numConfAtPos[res1];
    	int n2 = numConfAtPos[res2];
    	int firstIndex = pairwiseOffsets[getPairwiseIndex(res1, res2)];
    	for (int i1=0; i1<n1; i1++) {
    		for (int i2=0; i2<n2; i2++) {
    			pairwise.set(firstIndex + i1*n2 + i2, val.get(i1).get(i2));
    		}
    	}
    }
    
    public void setTupleValue(RCTuple tup, T val){
        //assign the given value to the specified RC tuple
        int tupSize = tup.pos.size();
        
        if(tupSize==1)//just a one-body quantity
            setOneBody( tup.pos.get(0), tup.RCs.get(0), val);
        else if(tupSize==2)//two-body
            setPairwise( tup.pos.get(0), tup.RCs.get(0), tup.pos.get(1), tup.RCs.get(1), val );
        else if(tupSize>2){//higher-order
        	// TODO: lazily-allocate space for higher-order values
            setHigherOrder(tup,val);
        }
        else
            throw new UnsupportedOperationException( "ERROR: Not supporting tuple size " + tupSize );
    }
    
    
    public void setHigherOrder(RCTuple tup, T val){
        //set a higher-order term
        //we need all pairs contained in tup to know about it
        
        //loop over pairs
        for(int index1=0; index1<tup.pos.size(); index1++){
            for(int index2=0; index2<index1; index2++){
                
                int pos1 = tup.pos.get(index1);
                int rc1 = tup.RCs.get(index1);
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);
                
                //put tup into the HigherTupleFinder for this pair
                HigherTupleFinder<T> htf = getHigherOrderTerms(pos1,rc1,pos2,rc2);
                
                //create a HigherTupleFinder if there is none yet
                if(htf==null){
                    htf = new HigherTupleFinder<>(defaultHigherInteraction);
                    setHigherOrderTerms(pos1, rc1, pos2, rc2, htf);
                }
                
                RCTuple subTup = tup.subtractMember(index1).subtractMember(index2);
                
                htf.setInteraction(subTup,val);
            }
        }
        
    }
    
    public HigherTupleFinder<T> getHigherOrderTerms(int res1, int conf1, int res2, int conf2) {
    	return higherTerms.get(getPairwiseIndex(res1, conf1, res2, conf2));
    }
    
    private void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<T> val) {
    	higherTerms.set(getPairwiseIndex(res1, conf1, res2, conf2), val);
    }
}
