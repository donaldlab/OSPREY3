package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;
import java.util.Iterator;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
public class PruningMatrixInverted extends PruningMatrix {

	private PruningMatrix other;
	private MSSearchProblem search;

	public PruningMatrixInverted(
			MSSearchProblem search,
			PruningMatrix other
			) {
		super();
		this.other = other;
		this.search = search;
	}

	private boolean somethingPrunedForAAType(int res, String rcAAType) {
		for(int index : other.prunedRCsAtPos(res)) {
			String rcAAType2 = search.confSpace.posFlex.get(res).RCs.get(index).AAType;
			if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
			else
				return true;
		}
		return false;
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		String rcAAType = search.confSpace.posFlex.get(res).RCs.get(index).AAType;

		//if not in specified list, then already marked as pruned in reduced matrix
		if(!search.settings.AATypeOptions.get(res).contains(rcAAType)) 
			return true;

		//if in specified aa list, invert.
		//except...if no rcs are pruned for a specified aa in our list, there 
		//will be no confs.
		else {
			if(somethingPrunedForAAType(res, rcAAType))
				return !other.getOneBody(res, index);

			//nothing pruned, so we need to keep this rc
			return false;
		}
	}

	private boolean contains(int res1, int index1, int res2, int index2) {
		String rcAAType1 = search.confSpace.posFlex.get(res1).RCs.get(index1).AAType;
		String rcAAType2 = search.confSpace.posFlex.get(res2).RCs.get(index2).AAType;
		boolean c1 = search.settings.AATypeOptions.get(res1).contains(rcAAType1);
		boolean c2 = search.settings.AATypeOptions.get(res2).contains(rcAAType2);
		return c1 && c2;
	}

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		if(contains(res1, index1, res2, index2)) {
			if(!getOneBody(res1, index1) && !getOneBody(res2, index2))
				return false;
		}
		return true;
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
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int conf1, int res2, int conf2) {
		throw new UnsupportedOperationException("higher order terms not yet supported with value inversion");
	}

	@Override
	public void setHigherOrderTerms(int res1, int conf1, int res2, int conf2, HigherTupleFinder<Boolean> val) {
		dontwrite();
	}

	@Override
	public int getNumPos() {
		return other.getNumPos();
	}

	@Override
	public int getNumConfAtPos(int pos) {
		return other.getNumConfAtPos(pos);
	}

	private void dontwrite() {
		throw new UnsupportedOperationException("this inverted pruning matrix is read-only");
	}

}
