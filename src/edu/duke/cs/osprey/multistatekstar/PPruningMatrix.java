package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Iterator;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Used for p* in full and partial sequences
 *
 */

@SuppressWarnings("serial")
public class PPruningMatrix extends QPruningMatrix {
	protected QPruningMatrix other;

	public PPruningMatrix(QPruningMatrix other){
		super(other);
		this.other = other;
		this.assignedAATypeOptions = other.assignedAATypeOptions;
		this.assignedFlexRes = other.assignedFlexRes;
		this.sp = other.sp;
		
		if(!isValid())
			throw new RuntimeException("ERROR: did not prune all RCs outside of assigned AA type options");
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		String rcAAType = other.sp.confSpace.posFlex.get(res).RCs.get(index).AAType;

		//if not in specified list, then already marked as pruned in reduced matrix
		if(!other.assignedAATypeOptions.get(res).contains(rcAAType)) 
			return true;

		//if in specified aa list, invert.
		//except...if no rcs are pruned for a specified aa in our list, there 
		//will be no P confs.
		else {
			if(somethingPrunedForAAType(res, rcAAType))
				return !other.getOneBody(res, index);

			//nothing pruned, so we need to keep this rc
			return false;
		}
	}

	public boolean somethingPrunedForAAType(int res, String rcAAType) {
		for(int index : other.prunedRCsAtPos(res)) {
			String rcAAType2 = other.sp.confSpace.posFlex.get(res).RCs.get(index).AAType;
			if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
			else
				return true;
		}
		return false;
	}

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		if(other.contains(res1, index1, res2, index2))
			return other.getPairwise(res1, index1, res2, index2);
		return true;
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		throw new UnsupportedOperationException("ERROR: higher order terms are not supported in P pruning matrix");
	}

	public PruningMatrix invert() {
		return new P2PruningMatrix(this);
	}

	public boolean isFullyDefined() {
		return other.isFullyDefined();
	}

	public PruningMatrix getParent() {
		return other.parent;
	}

	@Override
	public int getNumPos() {
		return other.getNumPos();
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
}
