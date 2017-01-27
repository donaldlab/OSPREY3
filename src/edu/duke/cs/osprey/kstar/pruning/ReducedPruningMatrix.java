package edu.duke.cs.osprey.kstar.pruning;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.kstar.KSSearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;

@SuppressWarnings("serial")
public class ReducedPruningMatrix extends PruningMatrix {

	// one body, pairwise, and higher order will be called by an object that sees
	// only the reduced number of positions. when accessing RCs in the pruneMat, 
	// it is necessary to convert from reduced position to the absolute position.
	// however, reducedAllowedAAs in sp has also been cut down, so we access that using 
	// the position specified by the caller

	protected UpdatedPruningMatrix upm;
	protected KSSearchProblem sp;

	public ReducedPruningMatrix(KSSearchProblem sp) {

		this.upm = new UpdatedPruningMatrix(sp.pruneMat);
		this.sp = sp;
		this.setPruningInterval(sp.pruneMat.getPruningInterval());
	}


	public ReducedPruningMatrix(KSSearchProblem sp, UpdatedPruningMatrix upm) {

		this.sp = sp;
		this.upm = upm;
		this.setPruningInterval(sp.pruneMat.getPruningInterval());
	}


	public UpdatedPruningMatrix getUpdatedPruningMatrix() {
		return upm;
	}


	@Override
	public Boolean getOneBody(int res, int index) {

		Integer pos = sp.posNums.get(res);

		return upm.getOneBody(pos, index);
	}


	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		return upm.getPairwise(pos1, index1, pos2, index2);
	}


	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		return upm.getHigherOrderTerms(pos1, index1, pos2, index2);
	}


	@Override
	public int getNumConfAtPos(int pos) {

		Integer pos1 = sp.posNums.get(pos);

		return sp.pruneMat.getNumConfAtPos(pos1);
	}


	@Override
	public int getNumPos() {
		return sp.posNums.size();
	}
	
	public int countUpdates() {
		return upm.countUpdates();
	}
}
