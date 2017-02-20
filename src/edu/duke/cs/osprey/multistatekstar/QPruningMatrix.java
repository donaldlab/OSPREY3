package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * Used for q*, q' in full and partial sequences
 *
 */
@SuppressWarnings("serial")
public class QPruningMatrix extends UpdatedPruningMatrix {

	protected int numPos;
	protected ArrayList<Integer> ignoredAbsPos;
	protected HashMap<Integer, Integer> red2AbsPos = new HashMap<>();
	protected HashMap<Integer, Integer> abs2RedPos = new HashMap<>();

	/**
	 * assuming parent is the original, all-seq matrix
	 * @param parentSp
	 * @param redFlexRes
	 * @param redAATypeOptions
	 */
	public QPruningMatrix(SearchProblem parentSp, 
			ArrayList<String> redFlexRes,
			ArrayList<ArrayList<String>> redAATypeOptions) {

		super(parentSp.pruneMat);
		numPos = redFlexRes.size();

		//create positions map
		int redPos = -1, absPos;//reduced pos, absolute pos
		for(String res : redFlexRes) {
			redPos++;
			absPos = parentSp.flexRes.indexOf(res);
			red2AbsPos.put(redPos, absPos);
			abs2RedPos.put(absPos, redPos);
		}

		ignoredAbsPos = new ArrayList<>();
		for(absPos=0;absPos<parent.getNumPos();++absPos) {
			if(!abs2RedPos.containsKey(absPos))
				ignoredAbsPos.add(absPos);
		}
		ignoredAbsPos.trimToSize();
		
		markNonAATypeOptionsAsUnPruned(parentSp, redAATypeOptions);
	}

	/**
	 * mark as pruned those rcs not corresponding to the reduced aa options
	 * @param parentSp
	 * @param redAATypeOptions
	 * @return
	 */
	private void markNonAATypeOptionsAsUnPruned(SearchProblem parentSp, 
			ArrayList<ArrayList<String>> redAATypeOptions) {
		//assuming parent is the original pruning matrix!
		Integer redPos;//reducedpos
		for(int absPos=0;absPos<parent.getNumPos();++absPos) {
			for(int rc : parent.unprunedRCsAtPos(absPos)) {
				String rcAAType = parentSp.confSpace.posFlex.get(absPos).RCs.get(rc).AAType;
				//not in reduced position, not a desired AA type
				if((redPos = abs2RedPos(absPos))==null || !redAATypeOptions.get(redPos).contains(rcAAType))
					super.markAsPruned(new RCTuple(absPos, rc));
			}
		}

		if(!isValid(parentSp, redAATypeOptions))
			throw new RuntimeException("ERROR: did not prune all RCs outside of reduced AA type options");
	}

	private boolean isValid(SearchProblem parentSp, 
			ArrayList<ArrayList<String>> redAATypeOptions) {
		//iterate through all redPos to get unpruned AAs
		ArrayList<ArrayList<String>> aaTypeOptions = new ArrayList<>();
		for(int redPos=0;redPos<numPos;++redPos) {
			aaTypeOptions.add(new ArrayList<>());
			for(int rc : unprunedRCsAtPos(redPos)) {
				String rcAAType = parentSp.confSpace.posFlex.get(red2AbsPos(redPos)).RCs.get(rc).AAType;
				if(!aaTypeOptions.get(redPos).contains(rcAAType))
					aaTypeOptions.get(redPos).add(rcAAType);
			}
		}

		for(int redPos=0;redPos<numPos;++redPos) {
			ArrayList<String> aasAtRedPos = new ArrayList<>(redAATypeOptions.get(redPos));
			Collections.sort(aasAtRedPos);
			Collections.sort(aaTypeOptions.get(redPos));

			if(!aasAtRedPos.equals(aaTypeOptions.get(redPos)))
				return false;
		}

		return true;
	}
	
	public boolean isFullyDefined() {
		return ignoredAbsPos.size()==0;
	}

	public ArrayList<Integer> getIgnoredAbsPos() {
		return ignoredAbsPos;
	}

	public PruningMatrix getParent() {
		return parent;
	}

	@Override
	public int getNumPos() {
		return numPos;
	}

	protected Integer red2AbsPos(int redPos) {
		return red2AbsPos.get(redPos);//can be null
		//Integer ans = red2AbsPos.get(redPos);
		//if(ans==null) throw new RuntimeException("ERROR: there is no position mapping for "+redPos);
		//return ans;
	}

	protected Integer abs2RedPos(int absPos) {
		return abs2RedPos.get(absPos);//can be null
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		return super.getOneBody(red2AbsPos(res), index);
	}

	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		return super.getPairwise(red2AbsPos(res1), index1, red2AbsPos(res2), index2);
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		return super.getHigherOrderTerms(red2AbsPos(res1), index1, red2AbsPos(res2), index2);
	}

	@Override
	public int getNumConfAtPos(int pos) {
		return super.getNumConfAtPos(red2AbsPos(pos));
	}

	public BigInteger getNumReducedUnprunedConfs() {
		BigInteger ans = BigInteger.ONE;
		for(int redPos=0;redPos<numPos;++redPos) {
			ans = ans.multiply(BigInteger.valueOf(unprunedRCsAtPos(redPos).size()));
			if(ans.compareTo(BigInteger.ZERO)==0) 
				return ans;
		}
		return ans;
	}
	
	/**
	 * returns unpruned rcs at all positions: reduced rcs at relevant positions
	 * and parent unpruned rcs at other positions
	 * @return
	 */
	public ArrayList<ArrayList<Integer>> getAllUnprunedRCsByAbsPos() {
		ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
		for(int absPos=0;absPos<parent.getNumPos();++absPos) {
			ans.add(new ArrayList<>());
			Integer redPos;
			if((redPos=abs2RedPos.get(absPos))==null)
				//this position is not in the reduced set, so get parent RCs
				ans.get(absPos).addAll(parent.unprunedRCsAtPos(absPos));
			else
				//get reduced rcs
				ans.get(absPos).addAll(super.unprunedRCsAtPos(redPos));
		}
		return ans;
	}
}
