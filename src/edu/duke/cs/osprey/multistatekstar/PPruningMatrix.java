package edu.duke.cs.osprey.multistatekstar;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Used for p* in full and partial sequences
 *
 */

@SuppressWarnings("serial")
public class PPruningMatrix extends QPruningMatrix {
	QPruningMatrix other;
	
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

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		Integer absPos1 = other.red2AbsPos(res1), absPos2 = other.red2AbsPos(res2);
		String rcAAType1 = other.sp.confSpace.posFlex.get(absPos1).RCs.get(index1).AAType;
		String rcAAType2 = other.sp.confSpace.posFlex.get(absPos2).RCs.get(index2).AAType;

		//if in specified aa list, invert
		if(other.redAATypeOptions.get(res1).contains(rcAAType1) && other.redAATypeOptions.get(res2).contains(rcAAType2))
			return !other.getPairwise(res1, index1, res2, index2);

		//not in sequence(s) of interest. must be already pruned in Q matrix
		return true;
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		throw new UnsupportedOperationException("ERROR: higher order terms not supported in P pruning matrix");
	}
	
	public PruningMatrix invert() {
		return other;
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
}
