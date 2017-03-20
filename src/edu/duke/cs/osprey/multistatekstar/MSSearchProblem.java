package edu.duke.cs.osprey.multistatekstar;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningMethod;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 * 
 */
@SuppressWarnings("serial")
public class MSSearchProblem extends SearchProblem {

	MSSearchSettings settings;
	private int numDefinedPos;

	public MSSearchProblem(SearchProblem other, 
			MSSearchSettings settings) {
		super(other);
		if(settings==null) throw new RuntimeException("ERROR: search settings cannot be null");
		this.settings = settings;
		this.pruneMat = getReducedPruningMatrix();
		this.allowedAAs = settings.AATypeOptions;
		this.flexRes = settings.mutRes;//-1 for unassigned positions
		this.numDefinedPos = other.confSpace.numPos-Collections.frequency(flexRes, "-1");
	}
	
	public ArrayList<Integer> getPos(boolean defined) {
		ArrayList<Integer> ans = new ArrayList<>();
		for(int i=0;i<flexRes.size();++i) {
			if(!defined && flexRes.get(i).equals("-1")) ans.add(i);//get undefined pos
			else if(defined && !flexRes.get(i).equals("-1")) ans.add(i);//get defined pos
		}
		ans.trimToSize();
		return ans;
	}
	
	public int getNumDefinedPos() {
		return numDefinedPos;
	}
	
	public boolean isFullyDefined() {
		return numDefinedPos==confSpace.numPos;
	}
	
	public ArrayList<Integer> unprunedAtPos(QPruningMatrix pruneMat, int pos, String AAType) {
		ArrayList<Integer> ans = new ArrayList<>();
		for(int rc : pruneMat.unprunedRCsAtPos(pos)) {
			String type = confSpace.posFlex.get(pos).RCs.get(rc).AAType;
			if(!AAType.equalsIgnoreCase(type)) continue;
			ans.add(rc);
		}
		return ans;
	}

	public QPruningMatrix prunePmat() {
		return prunePmat(this, settings.stericThreshold, settings.stericThreshold);
	}
	
	private QPruningMatrix getReducedPruningMatrix() {
		return new QPruningMatrix(this, settings.mutRes, settings.AATypeOptions);
	}
	
	private QPruningMatrix prunePmat(SearchProblem search, double pruningWindow, double stericThresh) {

		QPruningMatrix qpm = (QPruningMatrix) pruneMat;
		//don't want to overprune
		BigInteger minUnprunedConfs = BigInteger.valueOf(65536);
		
		//single sequence type dependent pruning for better efficiency
		//now do any consequent singles & pairs pruning
		int numUpdates = qpm.countUpdates();
		int oldNumUpdates;

		Pruner dee = new Pruner(search, qpm, true, stericThresh, pruningWindow, 
				search.useEPIC, search.useTupExpForSearch);
		dee.setVerbose(false);

		do {//repeat as long as we're pruning things
			oldNumUpdates = numUpdates;
			dee.prune("GOLDSTEIN");
			//pairs pruning can take a LONG time
			if(!contSCFlex || 
					(contSCFlex && dee.enumerateCandidates(PruningMethod.getMethod("GOLDSTEIN PAIRS FULL")).size() < 32768))
				dee.prune("GOLDSTEIN PAIRS FULL");
			numUpdates = qpm.countUpdates();
			
		} while (numUpdates > oldNumUpdates && 
				qpm.getNumUnprunedConfs().compareTo(minUnprunedConfs) > 0);
		
		return qpm;
	}
}
