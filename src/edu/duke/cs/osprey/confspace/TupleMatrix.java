/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.function.BiConsumer;

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


    // default implementation of tuple trees for higher order tuples

    default boolean hasHigherOrderTuples() {
    	return false;
	}

	/**
	 * returns the TupleTree containing the higher order tuples, or null if none
	 *
	 * use this for reading tuples from the matrix
	 */
	default TupleTree<T> getHigherOrderTuples(int pos1, int rc1, int pos2, int rc2) {
		throw new UnsupportedOperationException();
	}

	/**
	 * returns the TupleTree containing the higher order tuples
	 *
	 * if a tuple tree doesn't exist for this pair yet, a new one is created
	 *
	 * the tree can be used to put/get individual tuples
	 *
	 * use this for writing tuples to the matrix
	 */
	default TupleTree<T> getOrMakeHigherOrderTuples(int pos1, int rc1, int pos2, int rc2) {
		throw new UnsupportedOperationException();
	}

	default T getTuple(RCTuple tuple) {
		switch (tuple.size()) {

			case 1:
				return getOneBody(
					tuple.pos.get(0), tuple.RCs.get(0)
				);

			case 2:
				return getPairwise(
					tuple.pos.get(0), tuple.RCs.get(0),
					tuple.pos.get(1), tuple.RCs.get(1)
				);

			default:
				TupleTree<T> tree = getHigherOrderTuples(
					tuple.pos.get(0), tuple.RCs.get(0),
					tuple.pos.get(1), tuple.RCs.get(1)
				);
				if (tree != null) {
					return tree.get(tuple);
				}
				return null;
		}
	}

	default void setTuple(RCTuple tuple, T val) {
		switch (tuple.size()) {

			case 1:
				setOneBody(
					tuple.pos.get(0), tuple.RCs.get(0),
					val
				);
			break;

			case 2:
				setPairwise(
					tuple.pos.get(0), tuple.RCs.get(0),
					tuple.pos.get(1), tuple.RCs.get(1),
					val
				);
			break;

			default:
				getOrMakeHigherOrderTuples(
					tuple.pos.get(0), tuple.RCs.get(0),
					tuple.pos.get(1), tuple.RCs.get(1)
				).put(tuple, val);
			break;
		}
	}

	/**
	 * iterate over all higher-order (n>2) tuples matching the conformation
	 */
	default void forEachHigherOrderTupleIn(int[] conf, BiConsumer<RCTuple,T> callback) {

		if (hasHigherOrderTuples()) {

			int numPos = getNumPos();
			for (int pos1=2; pos1<numPos; pos1++) {

				int rc1 = conf[pos1];
				if (rc1 == Conf.Unassigned) {
					continue;
				}

				for (int pos2=1; pos2<pos1; pos2++) {

					int rc2 = conf[pos2];
					if (rc2 == Conf.Unassigned) {
						continue;
					}

					TupleTree<T> tree = getHigherOrderTuples(pos1, rc1, pos2, rc2);
					if (tree != null) {
						tree.forEachIn(conf, callback);
					}
				}
			}
		}
	}
}
