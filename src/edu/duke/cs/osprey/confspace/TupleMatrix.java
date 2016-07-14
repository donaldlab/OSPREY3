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
    
    double getPruningInterval();
    
    int getNumPos();
    int getNumConfAtPos(int pos);
    
    T getOneBody(int res, int conf);
    void setOneBody(int res, int conf, T val);
    void setOneBody(int res, ArrayList<T> val);
    
    T getPairwise(int res1, int conf1, int res2, int conf2);
    void setPairwise(int res1, int conf1, int res2, int conf2, T val);
    void setPairwise(int res1, int res2, ArrayList<ArrayList<T>> val);
    
    boolean hasHigherOrderTerms();
    void setTupleValue(RCTuple tup, T val);
    HigherTupleFinder<T> getHigherOrderTerms(int res1, int conf1, int res2, int conf2);
    void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<T> val);
}
