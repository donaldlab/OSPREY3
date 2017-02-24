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
		if(!isValid(other.sp, other.redAATypeOptions))
			throw new RuntimeException("ERROR: did not prune all RCs outside of reduced AA type options");
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		Integer absPos = other.red2AbsPos(res);
		String rcAAType = other.sp.confSpace.posFlex.get(absPos).RCs.get(index).AAType;

		//if not in specified list, then already marked as pruned in reduced matrix
		if(!other.redAATypeOptions.get(res).contains(rcAAType)) 
			return true;

		//if in specified aa list, invert.
		//except...if no rcs are pruned for a specified aa in our list, there 
		//will be no P confs.
		else {
			int numPrunedForAAType = 0;//count all pruned aas of type res-index
			for(int index2 : other.prunedRCsAtPos(res)) {
				String rcAAType2 = other.sp.confSpace.posFlex.get(absPos).RCs.get(index2).AAType;
				if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
				else {
					numPrunedForAAType++;
					break;
				}
			}

			if(numPrunedForAAType == 0)
				return other.getOneBody(res, index);

			return !other.getOneBody(res, index);
		}
	}

	public boolean somethingPrunedForAAType(String rcAAType, int res) {
		Integer absPos = other.red2AbsPos(res);
		for(int index : other.prunedRCsAtPos(res)) {
			String rcAAType2 = other.sp.confSpace.posFlex.get(absPos).RCs.get(index).AAType;
			if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
			else
				return true;
		}
		return false;
	}

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		//AAO: this logic is very tricky...this call is always preceded by a call to
		//either (un)prunedrcsatpos(res1) and (un)prunedrcsatpos(res2), which 
		//does the inversion!, so the correct thing to do is to return getpairwise
		//...also, the call to (un)prunedrcsatpos will have restricted the rcs to
		//the AAs of interest.
		return other.getPairwise(res1, index1, res2, index2);
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		throw new UnsupportedOperationException("ERROR: higher order terms are not supported in P pruning matrix");
	}

	public PruningMatrix invert() {
		return new P2PruningMatrix(this);
	}

	public boolean isFullyDefined() {
		return other.ignoredAbsPos.size()==0;
	}

	protected Integer red2AbsPos(int redPos) {
		return other.red2AbsPos.get(redPos);//can be null
	}

	protected Integer abs2RedPos(int absPos) {
		return other.abs2RedPos.get(absPos);//can be null
	}

	public ArrayList<Integer> getIgnoredAbsPos() {
		return other.ignoredAbsPos;
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
