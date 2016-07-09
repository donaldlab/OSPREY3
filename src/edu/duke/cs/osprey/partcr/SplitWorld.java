package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCIndexMap;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.LazyEnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.ShellDistribution;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class SplitWorld {
	
	private SearchProblem search;
	private RCSplits splits;
	private List<RCIndexMap> rcMaps;
	private SimpleEnergyCalculator ecalc;
	
	public SplitWorld(SearchProblem search, EnergyFunctionGenerator egen, ShellDistribution dist) {
		
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
		
		this.ecalc = new SimpleEnergyCalculator(egen, this.search.confSpace, search.shellResidues, dist);
		this.search.fullConfE = egen.fullConfEnergy(this.search.confSpace, search.shellResidues);
		this.search.emat = new LazyEnergyMatrix(search.emat, this.ecalc);
		this.search.pruneMat = new PruningMatrix(search.pruneMat);
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
	
	public void updateMatrices(double minBoundEnergy, double bestMinimizedEnergy, double Ew) {
	
		int numPos = search.confSpace.numPos;
		
		assert (search.emat != null);
		assert (search.pruneMat != null);
		
		// copy as many energies as possible from the old matrix
		// to a new matrix that is sized for the new conf space (with the added rcs)
		// for the rest of the energies, don't calculate them now. Instead, rely on lazy energy calculation
		// when combined with DEE pruning, this should cut down on the number of energies we have to calculate
		// and should make PartCR significantly faster
		
		double pruningInterval = bestMinimizedEnergy - minBoundEnergy + Ew;
		LazyEnergyMatrix oldEmat = (LazyEnergyMatrix)search.emat;
		LazyEnergyMatrix newEmat = new LazyEnergyMatrix(search.confSpace, pruningInterval, ecalc);
		PruningMatrix oldPmat = search.pruneMat;
		PruningMatrix newPmat = new PruningMatrix(search.confSpace, pruningInterval);
		
		for (int pos1=0; pos1<numPos; pos1++) {
			RCIndexMap rcMap1 = rcMaps.get(pos1);
			for (int rc1=0; rc1<newEmat.getNumConfAtPos(pos1); rc1++) {

				// get the old rc1, if any
				Integer oldRc1 = rc1;
				if (rcMap1 != null) {
					oldRc1 = rcMap1.newToOld(rc1);
				}
				
				// one-body
				if (oldRc1 != null) {
					
					// copy values from the old matrices
					if (oldEmat.hasOneBody(pos1, oldRc1)) {
						newEmat.setOneBody(pos1, rc1, oldEmat.getOneBody(pos1, oldRc1));
					}
					newPmat.setOneBody(pos1, rc1, oldPmat.getOneBody(pos1, oldRc1));
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
						
						if (oldRc1 != null && oldRc2 != null) {
							
							// copy values from the old matrices
							if (oldEmat.hasPairwise(pos1, oldRc1, pos2, oldRc2)) {
								newEmat.setPairwise(pos1, rc1, pos2, rc2, oldEmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
							}
							newPmat.setPairwise(pos1, rc1, pos2, rc2, oldPmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
						}
					}
				}
			}
		}
		
		search.emat = newEmat;
		search.pruneMat = newPmat;
		
		// do DEE to update the pruning matrix
		boolean typeDep = false;
		double boundsThresh = 100;
		int algOption = 1;
		boolean useFlags = true;
		boolean useTriples = false;
		boolean preDACS = false;
		boolean useEpic = false;
		boolean useTupExp = false;
		double stericThresh = 100;
		new PruningControl(
			search, pruningInterval, typeDep, boundsThresh, algOption,
			useFlags, useTriples, preDACS, useEpic, useTupExp, stericThresh
		).prune(PruningControl.ReportMode.None);
		
		// clear the maps now that we're done with them
		for (int pos1=0; pos1<numPos; pos1++) {
			rcMaps.set(pos1, null);
		}
	}
	
	public RCs makeRCs(int[] conf) {
		
		// make a pruning matrix that only leaves the parent rc at unsplit positions
		// and the split rcs at the split positions
		PruningMatrix pruneMat = new PruningMatrix(search.confSpace, 0);
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			
			List<RC> rcsAtPos = search.confSpace.posFlex.get(pos).RCs;
			RCSplits.RCInfo info = splits.getRCInfo(pos, conf[pos]);
			
			if (info.isSplit()) {
				
				// prune all but the splits
				for (int rc=0; rc<rcsAtPos.size(); rc++) {
					if (!info.isChild(rc)) {
						pruneMat.setOneBody(pos, rc, true);
					}
				}
				
			} else {
				
				// prune all but the parent
				for (int rc=0; rc<rcsAtPos.size(); rc++) {
					if (!info.isParent(rc)) {
						pruneMat.setOneBody(pos, rc, true);
					}
				}
			}
		}
		
		return new RCs(pruneMat);
	}

	public double improveBound(int[] conf) {
		RCs subRcs = makeRCs(conf);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, subRcs);
		return tree.nextLeafNode().getGScore();
	}
}
