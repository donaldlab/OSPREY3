package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
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

	public SearchProblem sp;
	public ArrayList<String> redFlexRes;
	public ArrayList<ArrayList<String>> redAATypeOptions;

	private int numPos;
	protected ArrayList<Integer> ignoredAbsPos;
	protected HashMap<Integer, Integer> red2AbsPos;
	protected HashMap<Integer, Integer> abs2RedPos;

	/**
	 * assuming parent is the original, all-seq matrix
	 * @param sp
	 * @param redFlexRes
	 * @param redAATypeOptions
	 */
	public QPruningMatrix(SearchProblem sp, 
			ArrayList<String> redFlexRes,
			ArrayList<ArrayList<String>> redAATypeOptions) {

		super(sp.pruneMat);
		this.sp = sp;
		this.redFlexRes = redFlexRes;
		this.redAATypeOptions = redAATypeOptions;
		numPos = redFlexRes.size();

		//create positions map
		red2AbsPos = new HashMap<>();
		abs2RedPos = new HashMap<>();
		int redPos = -1, absPos;//reduced pos, absolute pos
		for(String res : redFlexRes) {
			redPos++;
			absPos = sp.flexRes.indexOf(res);
			red2AbsPos.put(redPos, absPos);
			abs2RedPos.put(absPos, redPos);
		}

		ignoredAbsPos = new ArrayList<>();
		for(absPos=0;absPos<parent.getNumPos();++absPos) {
			if(!abs2RedPos.containsKey(absPos))
				ignoredAbsPos.add(absPos);
		}
		ignoredAbsPos.trimToSize();

		markNonAATypeOptionsAsPruned(sp, redAATypeOptions);
	}

	public QPruningMatrix(QPruningMatrix other) {
		super(other.sp.pruneMat);
	}

	/**
	 * mark as pruned those rcs not corresponding to the reduced aa options
	 * @param sp
	 * @param redAATypeOptions
	 * @return
	 */
	private void markNonAATypeOptionsAsPruned(SearchProblem sp, 
			ArrayList<ArrayList<String>> redAATypeOptions) {
		//assuming parent is the original pruning matrix!
		Integer redPos;//reducedpos
		for(int absPos=0;absPos<parent.getNumPos();++absPos) {
			for(int rc : parent.unprunedRCsAtPos(absPos)) {
				String rcAAType = sp.confSpace.posFlex.get(absPos).RCs.get(rc).AAType;
				//not in reduced position, not a desired AA type
				if((redPos = abs2RedPos(absPos))==null || !redAATypeOptions.get(redPos).contains(rcAAType))
					super.markAsPruned(new RCTuple(absPos, rc));
			}
		}

		if(!isValid(sp, redAATypeOptions))
			throw new RuntimeException("ERROR: did not prune all RCs outside of reduced AA type options");
	}

	protected boolean isValid(SearchProblem sp, 
			ArrayList<ArrayList<String>> redAATypeOptions) {
		//iterate through all redPos to get unpruned AAs
		ArrayList<ArrayList<String>> aaTypeOptions = new ArrayList<>();
		for(int redPos=0;redPos<getNumPos();++redPos) {
			aaTypeOptions.add(new ArrayList<>());
			for(int rc : unprunedRCsAtPos(redPos)) {
				String rcAAType = sp.confSpace.posFlex.get(red2AbsPos(redPos)).RCs.get(rc).AAType;
				if(!aaTypeOptions.get(redPos).contains(rcAAType))
					aaTypeOptions.get(redPos).add(rcAAType);
			}
		}

		//since sometimes, all rotamers at an aa are pruned (i.e. steric pruning)
		//therefore, the condition we test for is that all confirmed reduced aas 
		//(if any) must be in desired aa list
		for(int redPos=0;redPos<getNumPos();++redPos) {
			ArrayList<String> aas = aaTypeOptions.get(redPos);
			for(String aa : aas) {
				if(!redAATypeOptions.get(redPos).contains(aa)) return false;
			}
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
		for(int redPos=0;redPos<getNumPos();++redPos) {
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
				ans.get(absPos).addAll(unprunedRCsAtPos(redPos));
		}
		return ans;
	}

	public PruningMatrix invert() {
		return new PPruningMatrix(this);
	}
}
