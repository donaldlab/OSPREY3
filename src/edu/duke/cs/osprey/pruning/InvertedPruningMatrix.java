package edu.duke.cs.osprey.pruning;

import java.util.ArrayList;
import java.util.Iterator;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;

public class InvertedPruningMatrix extends PruningMatrix {

	private static final long serialVersionUID = 1993173402264691539L;
	
	private PruningMatrix pmat;
	
	public InvertedPruningMatrix(PruningMatrix pmat) {
		this.pmat = pmat;
		
		// NOTE: don't call a super constructor so we don't allocate anything
		// we just want to interpose on calls to the other pmat
	}
	
	private Boolean invert(Boolean val) {
		if (val == null) {
			return null;
		}
		return !val;
	}
	
	@Override
	public Boolean getOneBody(int res, int conf) {
		return invert(pmat.getOneBody(res, conf));
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
	public Boolean getPairwise(int res1, int conf1, int res2, int conf2) {
		return invert(pmat.getPairwise(res1, conf1, res2, conf2));
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
	public void unprunedRCsAtPos(ArrayList<Integer> out, int pos) {
		pmat.unprunedRCsAtPos(out, pos);
	}
	
	@Override
	public ArrayList<Integer> unprunedRCsAtPos(int pos) {
		return pmat.unprunedRCsAtPos(pos);
	}
	
	@Override
	public void prunedRCsAtPos(ArrayList<Integer> out, int pos) {
		pmat.prunedRCsAtPos(out, pos);
	}
	
	@Override
	public ArrayList<Integer> prunedRCsAtPos(int pos) {
		return pmat.prunedRCsAtPos(pos);
	}
	
	@Override
	public ArrayList<RCTuple> unprunedRCTuplesAtPos(ArrayList<Integer> pos) {
		return pmat.unprunedRCTuplesAtPos(pos);
	}
	
	@Override
	public boolean isPruned(RCTuple tup) {
		return pmat.isPruned(tup);
	}
	
	@Override
	public boolean isPrunedHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Boolean> htf) {
		return pmat.isPrunedHigherOrder(tup, curIndex, htf);
	}
	
	@Override
	public void markAsPruned(RCTuple tup) {
		dontwrite();
	}
	
	@Override
	public int countPrunedRCs() {
		return pmat.countPrunedRCs();
	}
	
	@Override
	public int countPrunedPairs() {
		return pmat.countPrunedPairs();
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
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int conf1, int res2, int conf2) {
		throw new UnsupportedOperationException("higher order terms not yet supported with value inversion");
	}
	
	@Override
	public void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<Boolean> val) {
		dontwrite();
	}
	
	@Override
	public int getNumPos() {
		return pmat.getNumPos();
	}
	
	@Override
	public int getNumConfAtPos(int pos) {
		return pmat.getNumConfAtPos(pos);
	}
	
	private void dontwrite() {
		throw new UnsupportedOperationException("this inverted pruning matrix is read-only");
	}
}
