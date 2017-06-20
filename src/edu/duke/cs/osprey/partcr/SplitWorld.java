package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCIndexMap;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.LazyEnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class SplitWorld {
	
	private SearchProblem search;
	private RCSplits splits;
	private List<RCIndexMap> rcMaps;
	private SimpleEnergyCalculator ecalc;
	
	public SplitWorld(SearchProblem search, ForcefieldParams ffparams) {
		
		// copy the search problem, confspace, the positions, and the rcs
		this.search = new SearchProblem(search);
		this.search.confSpace = new ConfSpace(search.confSpace);
		for (int pos=0; pos<this.search.confSpace.numPos; pos++) {
			this.search.confSpace.posFlex.set(pos, new PositionConfSpace(search.confSpace.posFlex.get(pos)));
		}
		
		this.splits = new RCSplits(this.search.confSpace);
		
		this.rcMaps = new ArrayList<>();
		for (int pos=0; pos<this.search.confSpace.numPos; pos++) {
			this.rcMaps.add(null);
		}
		
		this.ecalc = new SimpleEnergyCalculator.Cpu(ffparams, this.search.confSpace, this.search.shellResidues);
		this.search.fullConfE = ecalc.getEnergyFunctionGenerator().fullConfEnergy(this.search.confSpace, this.search.shellResidues);
		this.search.emat = new LazyEnergyMatrix(this.search.emat, this.ecalc);
		this.search.pruneMat = new PruningMatrix(this.search.pruneMat);
	}
	
	public SearchProblem getSearchProblem() {
		return search;
	}
	
	public RCSplits getSplits() {
		return splits;
	}
	
	public RC getRC(int pos, int rc) {
		return search.confSpace.posFlex.get(pos).RCs.get(rc);
	}
	
	public void replaceRc(int pos, RC rc, List<RC> splitRCs) {
		
		// keep track of which children go to which parents
		splits.split(rc, splitRCs);
		
		// split the rc and save the index map for updateMatrices()
		RCIndexMap map = search.confSpace.posFlex.get(pos).replaceRC(rc, splitRCs);
		rcMaps.set(pos, map);
	}
	
	public RCIndexMap getRcMap(int pos) {
		return rcMaps.get(pos);
	}
	
	public void resizeMatrices() {
		
		int numPos = search.confSpace.numPos;
		
		assert (search.emat != null);
		
		// copy as many energies as possible from the old matrix
		// to a new matrix that is sized for the new conf space (with the added rcs)
		// for the rest of the energies, don't calculate them now. Instead, rely on lazy energy calculation
		// when combined with DEE pruning, this should cut down on the number of energies we have to calculate
		// and should make PartCR significantly faster
		
		LazyEnergyMatrix oldEmat = (LazyEnergyMatrix)search.emat;
		LazyEnergyMatrix newEmat = new LazyEnergyMatrix(search.confSpace, oldEmat.getPruningInterval(), ecalc);
		
		for (int pos1=0; pos1<numPos; pos1++) {
			RCIndexMap rcMap1 = rcMaps.get(pos1);
			for (int rc1=0; rc1<newEmat.getNumConfAtPos(pos1); rc1++) {

				// one-body
				
				// get the old rc1, if any
				Integer oldRc1 = rc1;
				if (rcMap1 != null) {
					oldRc1 = rcMap1.newToOld(rc1);
				}
				
				// are there old values to copy?
				if (oldRc1 != null) {
					
					// copy them
					if (oldEmat.hasOneBody(pos1, oldRc1)) {
						newEmat.setOneBody(pos1, rc1, oldEmat.getOneBody(pos1, oldRc1));
					}
					
				} else {
					
					// we'll definitely need the single energy, so just calculate it now
					newEmat.setOneBody(pos1, rc1, ecalc.calcSingle(pos1, rc1).energy);
				}
				
				// pairwise
				for (int pos2=0; pos2<pos1; pos2++) {
					RCIndexMap rcMap2 = rcMaps.get(pos2);
					for (int rc2=0; rc2<newEmat.getNumConfAtPos(pos2); rc2++) {
						
						// get the old rc2, if any
						Integer oldRc2 = rc2;
						if (rcMap2 != null) {
							oldRc2 = rcMap2.newToOld(rc2);
						}
						
						// are there old values to copy?
						if (oldRc1 != null && oldRc2 != null) {
							
							// copy them
							if (oldEmat.hasPairwise(pos1, oldRc1, pos2, oldRc2)) {
								newEmat.setPairwise(pos1, rc1, pos2, rc2, oldEmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
							}
							
						} else {
							
							// we might not need the pair energy, so don't calculate it yet
							// do it lazily later if we really need it
							// gives about a 2x speedup in practice
							//newEmat.setPairwise(pos1, rc1, pos2, rc2, ecalc.calcPair(pos1, rc1, pos2, rc2).getEnergy());
						}
					}
				}
			}
		}
		
		search.emat = newEmat;
		
		// clear the maps now that we're done with them
		for (int pos1=0; pos1<numPos; pos1++) {
			rcMaps.set(pos1, null);
		}
	}
	
	public ScoredConf translateConf(ScoredConf conf) {
		
		// do A* search to find the improved bound
		RCs subRcs = splits.makeRCs(conf.getAssignments());
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, subRcs)
			.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setNumIterations(1)
			).build();
		return tree.nextConf();
	}
}
