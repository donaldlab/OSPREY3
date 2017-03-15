package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Iterator;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Used as the pruned matrix in the second phase of K*.
 * Reports everything as pruned.
 *
 */
public class P2PruningMatrix extends PruningMatrix {

	PruningMatrix other;
	
	public P2PruningMatrix(PPruningMatrix other) {
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
	
	public PruningMatrix invert() {
		throw new UnsupportedOperationException("ERROR: inversion is not supported in P2 pruning matrix");
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
		other.unprunedRCsAtPos(out, pos);
	}
	
	@Override
	public ArrayList<Integer> unprunedRCsAtPos(int pos) {
		return other.unprunedRCsAtPos(pos);
	}
	
	@Override
	public void prunedRCsAtPos(ArrayList<Integer> out, int pos) {
		other.prunedRCsAtPos(out, pos);
	}
	
	@Override
	public ArrayList<Integer> prunedRCsAtPos(int pos) {
		return other.prunedRCsAtPos(pos);
	}
	
	@Override
	public ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos) {
		return other.unprunedRCTuplesAtPos(pos);
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
		return other.countPrunedRCs();
	}
	
	@Override
	public int countPrunedPairs() {
		return other.countPrunedPairs();
	}

}
