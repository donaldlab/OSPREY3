/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.function.BiConsumer;
import java.util.function.Consumer;

public class TupleMatrixBoolean extends AbstractTupleMatrix<Boolean> {
	
	private static final long serialVersionUID = -1286639255089978027L;
	
    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
	private BitSet oneBody; // indices: res1, RC1
	private BitSet pairwise; // indices: res1, res2, RC1, RC2 where res1>res2
	
	protected TupleMatrixBoolean() {
		// do nothing
		// apparently UpdatedPruningMatrix wants to override all the methods,
		// but not use any of the storage here
	}
	
	public TupleMatrixBoolean(TupleMatrixBoolean other) {
		super(other);
		this.oneBody = (BitSet)other.oneBody.clone();
		this.pairwise = (BitSet)other.pairwise.clone();
	}
    
    public TupleMatrixBoolean(ConfSpace cSpace, double pruningInterval, boolean defaultHigherInteraction) {
    	super(cSpace, pruningInterval, defaultHigherInteraction);
    }

    public TupleMatrixBoolean(SimpleConfSpace confSpace, double pruningInterval, boolean defaultHigherInteraction) {
		super(confSpace, pruningInterval, defaultHigherInteraction);
	}
    
    public TupleMatrixBoolean(int numPos, int[] numAllowedAtPos, double pruningInterval, boolean defaultHigherInteraction) {
    	super(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }
    
    @Override
    protected void allocate(int numOneBody, int numPairwise) {
        oneBody = new BitSet(numOneBody);
        pairwise = new BitSet(numPairwise);
    }
    
    @Override
    public Boolean getOneBody(int res, int conf) {
    	return oneBody.get(getOneBodyIndex(res, conf));
    }
    
    @Override
    public void setOneBody(int res, int conf, Boolean val) {
    	oneBody.set(getOneBodyIndex(res, conf), val);
    }
    
    @Override
    public void setOneBody(int res, ArrayList<Boolean> val) {
    	int n = getNumConfAtPos(res);
    	for (int i=0; i<n; i++) {
    		oneBody.set(getOneBodyIndex(res, i), val.get(i));
    	}
    }
    
    @Override
    public Boolean getPairwise(int res1, int conf1, int res2, int conf2) {
    	return pairwise.get(getPairwiseIndex(res1, conf1, res2, conf2));
    }
    
    @Override
    public void setPairwise(int res1, int conf1, int res2, int conf2, Boolean val) {
    	pairwise.set(getPairwiseIndex(res1, conf1, res2, conf2), val);
    }
    
    @Override
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<Boolean>> val) {
    	int n1 = getNumConfAtPos(res1);
    	int n2 = getNumConfAtPos(res2);
    	for (int i1=0; i1<n1; i1++) {
    		for (int i2=0; i2<n2; i2++) {
    			pairwise.set(getPairwiseIndex(res1, i1, res2, i2), val.get(i1).get(i2));
    		}
    	}
    }

	@Override
	public String toString() {
		return toString(6, (isPruned) -> {
			return isPruned ? "  X   " : "      ";
		});
	}
}
