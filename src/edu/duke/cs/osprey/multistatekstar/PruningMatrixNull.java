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

package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Iterator;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * returns everything as pruned
 *
 */

@SuppressWarnings("serial")
public class PruningMatrixNull extends PruningMatrix {

	private PruningMatrix other;
	
	public PruningMatrixNull(PruningMatrix other) {
		super();
		this.other = other;
	}
	
	@Override
	public Boolean getOneBody(int res, int index) {
		return true;
	}
	
	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		return true;
	}
	
	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		throw new UnsupportedOperationException("ERROR: higher order terms are not supported in P2 pruning matrix");
	}
	
	@Override
	public int getNumPos() {
		return other.getNumPos();
	}
	
	@Override
	public int getNumConfAtPos(int pos) {
		return 0;
	}
	
	@Override
	public void setOneBody(int res, int conf, Boolean val) {
		dontwrite();
	}

	@Override
	public void setOneBody(int res, ArrayList<Boolean> val) {
		dontwrite();
	}

	@Override
	public void setPairwise(int res1, int conf1, int res2, int conf2, Boolean val) {
		dontwrite();
	}

	@Override
	public void setPairwise(int res1, int res2, ArrayList<ArrayList<Boolean>> val) {
		dontwrite();
	}

	@Override
	public void markAsPruned(RCTuple tup) {
		dontwrite();
	}

	@Override
	public void fill(Boolean val) {
		dontwrite();
	}

	@Override
	public void fill(Iterator<Boolean> val) {
		dontwrite();
	}

	@Override
	public void setTupleValue(RCTuple tup, Boolean val) {
		dontwrite();
	}

	@Override
	public void setHigherOrder(RCTuple tup, Boolean val) {
		dontwrite();
	}

	@Override
	public void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<Boolean> val) {
		dontwrite();
	}

	private void dontwrite() {
		throw new UnsupportedOperationException("ERROR: P pruning matrix is read-only");
	}
	
	@Override
	public void unprunedRCsAtPos(ArrayList<Integer> out, int pos) {
		return;
	}
	
	@Override
	public ArrayList<Integer> unprunedRCsAtPos(int pos) {
		return new ArrayList<Integer>();
	}
	
	@Override
	public void prunedRCsAtPos(ArrayList<Integer> out, int pos) {
		return;
	}
	
	@Override
	public ArrayList<Integer> prunedRCsAtPos(int pos) {
		return new ArrayList<Integer>();
	}
	
	@Override
	public ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos) {
		return new ArrayList<RCTuple>();
	}
	
	@Override
	public boolean isPruned(RCTuple tup) {
		return true;
	}
	
	@Override
	public boolean isPrunedHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Boolean> htf) {
		return true;
	}
	
	@Override
	public int countPrunedRCs() {
		return 0;
	}
	
	@Override
	public int countPrunedPairs() {
		return 0;
	}
}
