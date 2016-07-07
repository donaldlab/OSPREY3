package edu.duke.cs.osprey.partcr;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
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
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.ShellDistribution;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public class SplitWorld {
	
	private ConfSpace confSpace;
	private RCSplits splits;
	private List<RCIndexMap> rcMaps;
	private EnergyMatrix emat;
	private SimpleEnergyCalculator ecalc;
	private EnergyFunction efunc;
	
	public SplitWorld(SearchProblem search, EnergyFunctionGenerator egen, ShellDistribution dist) {
		
		// copy the confspace, the positions, and the rcs
		this.confSpace = new ConfSpace(search.confSpace);
		for (int pos=0; pos<this.confSpace.numPos; pos++) {
			this.confSpace.posFlex.set(pos, new PositionConfSpace(search.confSpace.posFlex.get(pos)));
		}
		
		this.splits = new RCSplits(confSpace);
		
		this.rcMaps = new ArrayList<>();
		for (int pos=0; pos<this.confSpace.numPos; pos++) {
			this.rcMaps.add(null);
		}
		
		this.emat = new EnergyMatrix(search.emat);
		this.ecalc = new SimpleEnergyCalculator(egen, confSpace, search.shellResidues, dist);
		this.efunc = egen.fullConfEnergy(confSpace, search.shellResidues);
	}
	
	public ConfSpace getConfSpace() {
		return confSpace;
	}
	
	public EnergyMatrix getEmat() {
		return emat;
	}
	
	public RCSplits getSplits() {
		return splits;
	}
	
	public RC getRC(int pos, int rc) {
		return confSpace.posFlex.get(pos).RCs.get(rc);
	}
	
	public void replaceRc(int pos, RC rc, List<RC> splitRCs) {
		
		// keep track of which children go to which parents
		splits.split(rc, splitRCs);
		
		// split the rc and save the index map for updateMatrices()
		RCIndexMap map = confSpace.posFlex.get(pos).replaceRC(rc, splitRCs);
		rcMaps.set(pos, map);
	}
	
	public RCIndexMap getRcMap(int pos) {
		return rcMaps.get(pos);
	}
	
	public void updateEnergyMatrix() {
	
		int numPos = confSpace.numPos;
		
		assert (emat != null);
		
		// copy as many energies as possible from the old matrix
		// to a new matrix that is sized for the new conf space (with the added rotamers)
		EnergyMatrix newEmat = new EnergyMatrix(confSpace, emat.getPruningInterval());
		
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
					newEmat.setOneBody(pos1, rc1, emat.getOneBody(pos1, oldRc1));
					
				} else {
					
					// calculate new values
					newEmat.setOneBody(pos1, rc1, ecalc.calcSingle(pos1, rc1).getEnergy());
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
							newEmat.setPairwise(pos1, rc1, pos2, rc2, emat.getPairwise(pos1, oldRc1, pos2, oldRc2));
						
						} else {
							
							// calculate new values
							newEmat.setPairwise(pos1, rc1, pos2, rc2, ecalc.calcPair(pos1, rc1, pos2, rc2).getEnergy());
						}
					}
				}
			}
		}
		
		// clear the maps now that we're done with them
		for (int pos1=0; pos1<numPos; pos1++) {
			rcMaps.set(pos1, null);
		}
			
		emat = newEmat;
	}
	
	public double minimize(int[] conf) {
		return confSpace.minimizeEnergy(conf, efunc, null);
	}
	
	public RCs makeRCs(ConfAStarNode leafNode) {
		
		int[] conf = new int[confSpace.numPos];
		leafNode.getConf(conf);
		
		// make a pruning matrix that only leaves the parent rc at unsplit positions
		// and the split rcs at the split positions
		PruningMatrix pruneMat = new PruningMatrix(confSpace, 0);
		for (int pos=0; pos<confSpace.numPos; pos++) {
			
			List<RC> rcsAtPos = confSpace.posFlex.get(pos).RCs;
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

	public ConfAStarNode improveNode(ConfAStarNode node) {
		RCs subRcs = makeRCs(node);
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, 1, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(emat), hscorer, subRcs);
		return tree.nextLeafNode();
	}
}
