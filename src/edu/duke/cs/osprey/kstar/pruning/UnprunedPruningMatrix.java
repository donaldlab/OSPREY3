package edu.duke.cs.osprey.kstar.pruning;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.kstar.KSSearchProblem;

@SuppressWarnings("serial")
public class UnprunedPruningMatrix extends ReducedPruningMatrix {

	// one body, pairwise, and higher order will be called by an object that sees
	// only the reduced number of positions. when accessing RCs in the pruneMat, 
	// it is necessary to convert from reduced position to the absolute position.
	// however, reducedAllowedAAs in sp has also been cut down, so we access that using 
	// the position specified by the caller

	public UnprunedPruningMatrix(KSSearchProblem sp, UpdatedPruningMatrix upm, double pruningInterval) {
		super(sp, upm);
		this.setPruningInterval(pruningInterval);
	}


	@Override
	public Boolean getOneBody(int res, int index) {

		Integer pos = sp.posNums.get(res);
		String rcAAType = sp.confSpace.posFlex.get(pos).RCs.get(index).AAType;

		// if in specified aa list, invert
		if(sp.reducedAllowedAAs.get(res).contains(rcAAType))
			return false;

		return true;
	}


	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		String rcAAType1 = sp.confSpace.posFlex.get(pos1).RCs.get(index1).AAType;
		String rcAAType2 = sp.confSpace.posFlex.get(pos2).RCs.get(index2).AAType;

		// if in specified aa list, invert
		if(sp.reducedAllowedAAs.get(res1).contains(rcAAType1) && sp.reducedAllowedAAs.get(res2).contains(rcAAType2))
			return false;

		return true;
	}


	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {

		Integer pos1 = sp.posNums.get(res1);
		Integer pos2 = sp.posNums.get(res2);

		String rcAAType1 = sp.confSpace.posFlex.get(pos1).RCs.get(index1).AAType;
		String rcAAType2 = sp.confSpace.posFlex.get(pos2).RCs.get(index2).AAType;

		// allowedAAs has been reduced, so use resX positions
		if(sp.reducedAllowedAAs.get(res1).contains(rcAAType1) && sp.reducedAllowedAAs.get(res2).contains(rcAAType2))
			return upm.getHigherOrderTerms(pos1, index1, pos2, index2);

		return null;
	}
}
