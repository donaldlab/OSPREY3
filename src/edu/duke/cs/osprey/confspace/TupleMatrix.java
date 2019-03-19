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

    default T getOneBody(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1) {
    	return getOneBody(pos1.index, rc1.index);
	}
	default void setOneBody(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, T val) {
		setOneBody(pos1.index, rc1.index, val);
	}

	T getPairwise(int res1, int conf1, int res2, int conf2);
    void setPairwise(int res1, int conf1, int res2, int conf2, T val);
    void setPairwise(int res1, int res2, ArrayList<ArrayList<T>> val);

	default T getPairwise(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2) {
		return getPairwise(pos1.index, rc1.index, pos2.index, rc2.index);
	}
	default void setPairwise(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2, T val) {
		setPairwise(pos1.index, rc1.index, pos2.index, rc2.index, val);
	}

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
	 *
	 * pos1 < pos2 must be the lowest-order positions in the tuple
	 */
	default TupleTree<T> getHigherOrderTuples(int pos1, int rc1, int pos2, int rc2) {
		throw new UnsupportedOperationException();
	}

	/**
	 * returns the TupleTree containing the higher order tuples
	 *
	 * if a tuple tree doesn't exist for this pair yet, a new one is created
	 *
	 * pos1 < pos2 must be the lowest-order positions in the tuple
	 */
	default TupleTree<T> getOrMakeHigherOrderTuples(int pos1, int rc1, int pos2, int rc2) {
		throw new UnsupportedOperationException();
	}

	default T getTuple(RCTuple tuple) {
		tuple.checkSortedPositions();
		switch (tuple.size()) {

			case 0: throw new IllegalArgumentException("zero-length tuple");

			case 1: {
				int pos1 = tuple.pos.get(0);
				int rc1 = tuple.RCs.get(0);
				return getOneBody(pos1, rc1);
			}

			case 2: {
				// choose pos1,pos2 such that pos1 < pos2
				int pos1 = tuple.pos.get(1);
				int rc1 = tuple.RCs.get(1);
				int pos2 = tuple.pos.get(0);
				int rc2 = tuple.RCs.get(0);
				return getPairwise(pos1, rc1, pos2, rc2);
			}

			default: {
				// choose pos1,pos2 such that pos1 < pos2 < pos3 ...
				int pos1 = tuple.pos.get(0);
				int rc1 = tuple.RCs.get(0);
				int pos2 = tuple.pos.get(1);
				int rc2 = tuple.RCs.get(1);
				TupleTree<T> tree = getHigherOrderTuples(pos1, rc1, pos2, rc2);
				if (tree != null) {
					return tree.get(tuple);
				}
				return null;
			}
		}
	}

	default void setTuple(RCTuple tuple, T val) {
		tuple.checkSortedPositions();
		switch (tuple.size()) {

			case 0: throw new IllegalArgumentException("zero-length tuple");

			case 1: {
				int pos1 = tuple.pos.get(0);
				int rc1 = tuple.RCs.get(0);
				setOneBody(pos1, rc1, val);
			} break;

			case 2: {
				// choose pos1,pos2 such that pos1 < pos2
				int pos1 = tuple.pos.get(1);
				int rc1 = tuple.RCs.get(1);
				int pos2 = tuple.pos.get(0);
				int rc2 = tuple.RCs.get(0);
				setPairwise(pos1, rc1, pos2, rc2, val);
			} break;

			default:
				// choose pos1,pos2 such that pos1 < pos2 < pos3 ...
				int pos1 = tuple.pos.get(0);
				int rc1 = tuple.RCs.get(0);
				int pos2 = tuple.pos.get(1);
				int rc2 = tuple.RCs.get(1);
				getOrMakeHigherOrderTuples(pos1, rc1, pos2, rc2).put(tuple, val);
			break;
		}
	}

	/**
	 * iterate over all higher-order (n>2) tuples matching the conformation
	 */
	default void forEachHigherOrderTupleIn(int[] conf, BiConsumer<RCTuple,T> callback) {

		if (!hasHigherOrderTuples()) {
			return;
		}

		int numPos = getNumPos();
		for (int pos1=1; pos1<numPos; pos1++) {

			int rc1 = conf[pos1];
			if (rc1 == Conf.Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

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

	/**
	 * iterate over all higher-order (n>2) tuples matching the conformation but containing the given position
	 */
	default void forEachHigherOrderTupleIn(int[] conf, int posa, BiConsumer<RCTuple,T> callback) {

		if (!hasHigherOrderTuples()) {
			return;
		}

		int numPos = getNumPos();
		for (int pos1=1; pos1<numPos; pos1++) {

			int rc1 = conf[pos1];
			if (rc1 == Conf.Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

				int rc2 = conf[pos2];
				if (rc2 == Conf.Unassigned) {
					continue;
				}

				TupleTree<T> tree = getHigherOrderTuples(pos1, rc1, pos2, rc2);
				if (tree != null) {
					tree.forEachIn(conf, posa, callback);
				}
			}
		}
	}

	/**
	 * iterate over all higher-order (n>2) tuples matching the conformation but containing the given positions
	 */
	default void forEachHigherOrderTupleIn(int[] conf, int posa, int posb, BiConsumer<RCTuple,T> callback) {

		if (!hasHigherOrderTuples()) {
			return;
		}


		int numPos = getNumPos();
		for (int pos1=1; pos1<numPos; pos1++) {

			int rc1 = conf[pos1];
			if (rc1 == Conf.Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

				int rc2 = conf[pos2];
				if (rc2 == Conf.Unassigned) {
					continue;
				}

				TupleTree<T> tree = getHigherOrderTuples(pos1, rc1, pos2, rc2);
				if (tree != null) {
					tree.forEachIn(conf, posa, posb, callback);
				}
			}
		}
	}
}
