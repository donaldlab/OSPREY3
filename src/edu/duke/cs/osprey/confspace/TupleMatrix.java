/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.Iterator;

/**
 *
 * @author mhall44
 */
public interface TupleMatrix<T> {
	
    void fill(T val);
    void fill(Iterator<T> val);
    
    public double getPruningInterval();
    
    public int getNumPos();
    public int getNumConfAtPos(int pos);
    
    public T getOneBody(int res, int conf);
    public void setOneBody(int res, int conf, T val);
    public void setOneBody(int res, ArrayList<T> val);
    
    public T getPairwise(int res1, int conf1, int res2, int conf2);
    public void setPairwise(int res1, int conf1, int res2, int conf2, T val);
    public void setPairwise(int res1, int res2, ArrayList<ArrayList<T>> val);
    
    public void setTupleValue(RCTuple tup, T val);
    public HigherTupleFinder<T> getHigherOrderTerms(int res1, int conf1, int res2, int conf2);
    public void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<T> val);
}
