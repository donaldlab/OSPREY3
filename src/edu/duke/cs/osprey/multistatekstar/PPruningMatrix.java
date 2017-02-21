package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Used for p* in full and partial sequences
 *
 */

@SuppressWarnings("serial")
public class PPruningMatrix extends QPruningMatrix {

	ArrayList<PositionConfSpace> posFlex;
	ArrayList<ArrayList<String>> redAATypeOptions;

	public PPruningMatrix(SearchProblem parentSp, ArrayList<String> redFlexRes,
			ArrayList<ArrayList<String>> redAATypeOptions) {
		super(parentSp, redFlexRes, redAATypeOptions);
		posFlex = parentSp.confSpace.posFlex;
		this.redAATypeOptions = redAATypeOptions;
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		Integer absPos = red2AbsPos(res);
		String rcAAType = posFlex.get(absPos).RCs.get(index).AAType;

		//if not in specified list, then already marked as pruned in reduced matrix
		if(!redAATypeOptions.get(res).contains(rcAAType)) 
			return true;

		//if in specified aa list, invert.
		//except...if no rcs are pruned for a specified aa in our list, there 
		//will be no P confs.
		else {
			int numPrunedForAAType = 0;//count all pruned aas of type res-index
			for(int index2 : prunedRCsAtPos(res)) {
				String rcAAType2 = posFlex.get(absPos).RCs.get(index2).AAType;
				if(!rcAAType2.equalsIgnoreCase(rcAAType)) continue;
				else {
					numPrunedForAAType++;
					break;
				}
			}

			if(numPrunedForAAType == 0)
				return super.getOneBody(res, index);

			return !super.getOneBody(res, index);
		}
	}

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		Integer absPos1 = red2AbsPos(res1), absPos2 = red2AbsPos(res2);
		String rcAAType1 = posFlex.get(absPos1).RCs.get(index1).AAType;
		String rcAAType2 = posFlex.get(absPos2).RCs.get(index2).AAType;

		//if in specified aa list, invert
		if(redAATypeOptions.get(res1).contains(rcAAType1) && redAATypeOptions.get(res2).contains(rcAAType2))
			return !super.getPairwise(res1, index1, res2, index2);

		//not in sequence(s) of interest. must be already pruned in Q matrix
		return true;
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		throw new UnsupportedOperationException("ERROR: higher order terms not supported in P pruning matrix");
	}
}
