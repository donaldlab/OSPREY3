package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;

import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
@SuppressWarnings("serial")
public class MSSearchProblem extends SearchProblem {

	public MSSearchSettings settings;
	private int numAssignedPos;

	public MSSearchProblem(SearchProblem other, 
			MSSearchSettings settings) {
		super(other);
		if(settings==null) throw new RuntimeException("ERROR: search settings cannot be null");
		this.settings = settings;
		//this.allowedAAs = settings.AATypeOptions;
		//this.flexRes = settings.mutRes;//-1 for unassigned positions
		this.numAssignedPos = other.confSpace.numPos-Collections.frequency(settings.mutRes, "-1");
	}

	public ArrayList<Integer> getPosNums(boolean assigned) {
		ArrayList<Integer> ans = new ArrayList<>();
		for(int i=0;i<settings.mutRes.size();++i) {
			if(!assigned && settings.mutRes.get(i).equals("-1")) ans.add(i);//get undefined pos
			else if(assigned && !settings.mutRes.get(i).equals("-1")) ans.add(i);//get defined pos
		}
		ans.trimToSize();
		return ans;
	}
	
	public ArrayList<String> getAssignedAAs() {
		ArrayList<String> ans = new ArrayList<>();
		for(int i : getPosNums(true)) {
			ans.add(settings.AATypeOptions.get(i).get(0));
		}
		ans.trimToSize();
		return ans;
	}

	public int getNumUnAssignedPos() {
		return confSpace.numPos-numAssignedPos;
	}
	
	public int getNumAssignedPos() {
		return numAssignedPos;
	}

	public boolean isFullyAssigned() {
		return numAssignedPos==confSpace.numPos;
	}

	public ArrayList<Integer> unprunedAtPos(PruningMatrix pruneMat, int pos, String AAType) {
		ArrayList<Integer> ans = new ArrayList<>();
		for(int rc : pruneMat.unprunedRCsAtPos(pos)) {
			String type = confSpace.posFlex.get(pos).RCs.get(rc).AAType;
			if(!AAType.equalsIgnoreCase(type)) continue;
			ans.add(rc);
		}
		return ans;
	}

	public void setPruningMatrix() {
		this.pruneMat = updatePruningMatrix(getPosNums(true), settings.AATypeOptions);
	}

	public void setPruningMatrix(PruningMatrix pmat) {
		this.pruneMat = pmat;
	}

	private BigInteger getNumConfs(PruningMatrix pmat) {
		BigInteger ans = BigInteger.ONE;
		for(int pos=0;pos<pmat.getNumPos();++pos) {
			ans = ans.multiply(BigInteger.valueOf(pmat.unprunedRCsAtPos(pos).size()));
			if(ans.compareTo(BigInteger.ZERO)==0) 
				return ans;
		}
		return ans;
	}

	public PruningMatrix updatePruningMatrix(
			ArrayList<Integer> splitPosNums, 
			ArrayList<ArrayList<String>> splitAAs
			) {
		UpdatedPruningMatrix ans = new UpdatedPruningMatrix(pruneMat);
		for(int pos : splitPosNums) {
			for(int rc : pruneMat.unprunedRCsAtPos(pos)) {
				String rcAAType = confSpace.posFlex.get(pos).RCs.get(rc).AAType;
				//not in reduced position, not a desired AA type
				if(!splitAAs.get(pos).contains(rcAAType))
					ans.markAsPruned(new RCTuple(pos, rc));
			}
		}
		return ans;
	}

	public PruningMatrix updatePruningMatrix(int splitPos, String splitAA) {
		UpdatedPruningMatrix ans = new UpdatedPruningMatrix(pruneMat);
		for(int rc : pruneMat.unprunedRCsAtPos(splitPos)) {
			String rcAAType = confSpace.posFlex.get(splitPos).RCs.get(rc).AAType;
			//not in reduced position, not a desired AA type
			if(!splitAA.equals(rcAAType))
				ans.markAsPruned(new RCTuple(splitPos, rc));
		}
		return ans;
	}

	private void prunePmat(SearchProblem search) {

		UpdatedPruningMatrix upmat = (UpdatedPruningMatrix) this.pruneMat;

		//don't want to overprune
		BigInteger minConfs = BigInteger.valueOf(65536);

		//single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = upmat.countUpdates();
		int oldNumUpdates;

		Pruner dee = new Pruner(search, upmat, true, settings.stericThreshold, 
				settings.pruningWindow, search.useEPIC, search.useTupExpForSearch);
		dee.setVerbose(false);

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = upmat.countUpdates();
		} while (numUpdates > oldNumUpdates && getNumConfs(upmat).compareTo(minConfs) > 0);
	}

	public void prunePmat() {
		prunePmat(this);
	}
}
