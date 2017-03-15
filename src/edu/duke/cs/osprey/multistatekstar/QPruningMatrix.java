package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
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
public class QPruningMatrix extends UpdatedPruningMatrix{

	public SearchProblem sp;
	public ArrayList<String> assignedFlexRes;
	public ArrayList<ArrayList<String>> assignedAATypeOptions;

	private int numPos;

	/**
	 * Create an update to the original pruning matrix to reflect the fact that
	 * only redflexres have assigned amino acid type(s). Unassigned residues in
	 * redflexres have a value of "-1". Corresponding values in redAAOptions are
	 * also "-1".
	 * @param sp
	 * @param flexRes
	 * @param AATypeOptions
	 */
	public QPruningMatrix(SearchProblem sp, 
			ArrayList<String> flexRes,
			ArrayList<ArrayList<String>> AATypeOptions) {

		super(sp.pruneMat);
		this.sp = sp;
		this.assignedFlexRes = flexRes;
		this.assignedAATypeOptions = AATypeOptions;
		numPos = flexRes.size();

		markNonAATypeOptionsAsPruned();
	}

	public QPruningMatrix(QPruningMatrix other) {
		super(other.sp.pruneMat);
	}

	/**
	 * mark as pruned those rcs not corresponding to the reduced aa options
	 * @return
	 */
	private void markNonAATypeOptionsAsPruned() {
		//assuming parent is the original pruning matrix!
		for(int pos=0;pos<parent.getNumPos();++pos) {
			for(int rc : parent.unprunedRCsAtPos(pos)) {
				String rcAAType = sp.confSpace.posFlex.get(pos).RCs.get(rc).AAType;
				//not in reduced position, not a desired AA type
				if(assignedFlexRes.get(pos).equals("-1") || !assignedAATypeOptions.get(pos).contains(rcAAType))
					super.markAsPruned(new RCTuple(pos, rc));
			}
		}

		if(!isValid())
			throw new RuntimeException("ERROR: did not prune all RCs outside of assigned AA type options");
	}

	protected boolean isValid() {
		//iterate through all redPos to get unpruned AAs
		ArrayList<ArrayList<String>> aaTypeOptions = new ArrayList<>();
		for(int pos=0;pos<getNumPos();++pos) {
			aaTypeOptions.add(new ArrayList<>());
			if(assignedFlexRes.get(pos).equals("-1")){
				aaTypeOptions.get(pos).addAll(assignedAATypeOptions.get(pos));
				continue;
			}
			for(int rc : unprunedRCsAtPos(pos)) {
				String rcAAType = sp.confSpace.posFlex.get(pos).RCs.get(rc).AAType;
				if(!aaTypeOptions.get(pos).contains(rcAAType))
					aaTypeOptions.get(pos).add(rcAAType);
			}
		}

		//sometimes, all rotamers at an aa are pruned (i.e. steric pruning)
		//therefore, the condition we test for is that all aas from reduced 
		//rotamer list (if any) must be in desired aa list
		for(int pos=0;pos<getNumPos();++pos) {
			ArrayList<String> aas = aaTypeOptions.get(pos);
			for(String aa : aas) {
				if(!assignedAATypeOptions.get(pos).contains(aa)) return false;
			}
		}

		return true;
	}
	
	public boolean contains(int res, int index) {
		String rcAAType = sp.confSpace.posFlex.get(res).RCs.get(index).AAType;
		boolean c = assignedAATypeOptions.get(res).contains(rcAAType);
		return c;
	}
	
	public boolean contains(int res1, int index1, int res2, int index2) {
		String rcAAType1 = sp.confSpace.posFlex.get(res1).RCs.get(index1).AAType;
		String rcAAType2 = sp.confSpace.posFlex.get(res2).RCs.get(index2).AAType;
		boolean c1 = assignedAATypeOptions.get(res1).contains(rcAAType1);
		boolean c2 = assignedAATypeOptions.get(res2).contains(rcAAType2);
		return c1 && c2;
	}

	public PruningMatrix getParent() {
		return parent;
	}

	@Override
	public int getNumPos() {
		return numPos;
	}

	@Override
	public Boolean getOneBody(int res, int index) {
		if(contains(res, index))
			return super.getOneBody(res, index);
		return true;
	}
	
	@Override
	public Boolean getPairwise(int res1, int index1, int res2, int index2) {
		if(contains(res1, index1, res2, index2))
			return super.getPairwise(res1, index1, res2, index2);
		return true;
	}

	@Override
	public HigherTupleFinder<Boolean> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
		if(contains(res1, index1, res2, index2))
			return super.getHigherOrderTerms(res1, index1, res2, index2);
		return null;
	}

	@Override
	public int getNumConfAtPos(int pos) {
		return super.getNumConfAtPos(pos);
	}

	public BigInteger getNumUnprunedConfs() {
		BigInteger ans = BigInteger.ONE;
		for(int pos=0;pos<getNumPos();++pos) {
			ans = ans.multiply(BigInteger.valueOf(unprunedRCsAtPos(pos).size()));
			if(ans.compareTo(BigInteger.ZERO)==0) 
				return ans;
		}
		return ans;
	}
	
	public int countUpdates() {
		return super.countUpdates();
	}

	public PruningMatrix invert() {
		return new PPruningMatrix(this);
	}
	
}
