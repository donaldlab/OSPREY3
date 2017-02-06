package edu.duke.cs.osprey.kstar.pruning;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.kstar.KSSearchProblem;

@SuppressWarnings("serial")
public class InvertedPruningMatrix extends ReducedPruningMatrix {

	// one body, pairwise, and higher order will be called by an object that sees
	// only the reduced number of positions. when accessing RCs in the pruneMat, 
	// it is necessary to convert from reduced position to the absolute position.
	// however, reducedAllowedAAs in sp has also been cut down, so we access that using 
	// the position specified by the caller

	public InvertedPruningMatrix(KSSearchProblem sp, UpdatedPruningMatrix upm) {
		super(sp, upm);
	}

	
	@Override
	public Boolean getOneBody(int res, int index) {

		Integer pos = sp.posNums.get(res);
		String rcAAType = sp.confSpace.posFlex.get(pos).RCs.get(index).AAType;

		// if in specified aa list, invert
		if(sp.reducedAllowedAAs.get(res).contains(rcAAType)) {

			int numPrunedForAAType = 0;
			for(int index2 : upm.prunedRCsAtPos(pos)) {
				String rcAAType2 = sp.confSpace.posFlex.get(pos).RCs.get(index2).AAType;
				if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
				else {
					numPrunedForAAType++;
					break;
				}
			}

			if(numPrunedForAAType == 0)
				return upm.getOneBody(pos,index);

			else
				return !upm.getOneBody(pos,index);
		}

		// not in sequence(s) of interest. must be already pruned in upm.
		// return true because we are always interested in unpruned confs at pos
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
			return !upm.getPairwise(pos1, index1, pos2, index2);

		// not in sequence(s) of interest. must be already pruned in upm.
		// return true because we are always interested in unpruned confs at pos
		return true;
	}


	// return null here for expediency
	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		//return parent.getHigherOrderTerms(res1, index1, res2, index2);
		return null;
	}
}
