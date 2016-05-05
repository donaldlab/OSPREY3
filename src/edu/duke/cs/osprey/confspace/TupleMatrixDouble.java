/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;

public class TupleMatrixDouble extends AbstractTupleMatrix<Double> {
	
	private static final long serialVersionUID = -1286639255089978027L;
	
    //note: tuples are sets not ordered pairs, i.e. E(i_r,j_s) = E(j_s,i_r), and pruning (i_r,j_s) means pruning (j_s,i_r)
	private double[] oneBody; // indices: res1, RC1
	private double[] pairwise; // indices: res1, res2, RC1, RC2 where res1>res2
    
    public TupleMatrixDouble(ConfSpace cSpace, double pruningInterval, double defaultHigherInteraction) {
    	super(cSpace, pruningInterval, defaultHigherInteraction);
    }
    
    public TupleMatrixDouble(int numPos, int[] numAllowedAtPos, double pruningInterval, double defaultHigherInteraction) {
    	super(numPos, numAllowedAtPos, pruningInterval, defaultHigherInteraction);
    }
    
    @Override
    protected void allocate(int numOneBody, int numPairwise) {
        oneBody = new double[numOneBody];
        pairwise = new double[numPairwise];
    }
    
    @Override
    public Double getOneBody(int res, int conf) {
    	return oneBody[getOneBodyIndex(res, conf)];
    }
    
    @Override
    public void setOneBody(int res, int conf, Double val) {
    	oneBody[getOneBodyIndex(res, conf)] = val;
    }
    
    @Override
    public void setOneBody(int res, ArrayList<Double> val) {
    	int n = getNumConfAtPos(res);
    	for (int i=0; i<n; i++) {
    		oneBody[getOneBodyIndex(res, i)] = val.get(i);
    	}
    }
    
    @Override
    public Double getPairwise(int res1, int conf1, int res2, int conf2) {
    	return pairwise[getPairwiseIndex(res1, conf1, res2, conf2)];
    }
    
    @Override
    public void setPairwise(int res1, int conf1, int res2, int conf2, Double val) {
    	pairwise[getPairwiseIndex(res1, conf1, res2, conf2)] = val;
    }
    
    @Override
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<Double>> val) {
    	int n1 = getNumConfAtPos(res1);
    	int n2 = getNumConfAtPos(res2);
    	for (int i1=0; i1<n1; i1++) {
    		for (int i2=0; i2<n2; i2++) {
    			pairwise[getPairwiseIndex(res1, i1, res2, i2)] = val.get(i1).get(i2);
    		}
    	}
    }
}
