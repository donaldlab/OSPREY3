package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

public interface ConfSearchFactory {
	
	ConfSearch make(EnergyMatrix emat, PruningMatrix pmat);
	
	// TODO: move this info some CFP-only place
	public static class Tools {
		
		public static ConfSearchFactory makeFromConfig(SearchProblem search, ConfigFileParser cfp) {
			return new ConfSearchFactory() {
				@Override
				public ConfSearch make(EnergyMatrix emat, PruningMatrix pmat) {
					
					if (emat.getNumPos() != pmat.getNumPos()) {
						throw new Error("energy matrix doesn't match pruning matrix, this is a bug");
					}
					
					if (search.searchNeedsHigherOrderTerms() || search.useEPIC) {
				
						// if we need higher-order or EPIC terms, use the old A* code
						return ConfTree.makeFull(search);
					}
					
					// when we don't need higher order terms, we can do fast pairwise-only things
					
					// get the appropriate configuration for A*
					AStarScorer gscorer = new PairwiseGScorer(emat);
					AStarScorer hscorer;
					AStarOrder order;
					RCs rcs = new RCs(pmat);
					
					// how many iterations of MPLP should we do?
					int numMPLPIters = cfp.getParams().getInt("NumMPLPIters");
					if (numMPLPIters <= 0) {
						
						// zero MPLP iterations is exactly the traditional heuristic, so use the fast implementation
						hscorer = new TraditionalPairwiseHScorer(emat, rcs);
						order = new DynamicHMeanAStarOrder();
						
					} else {
						
						// simple (but not exhaustive) testing showed node-based MPLP is basically
						// always faster than edge-based MPLP, so just use node-based MPLP all the time.
						// also, always use a static order with MPLP
						// MPLP isn't optimized to do differential node scoring quickly so DynamicHMean is super slow!
						double convergenceThreshold = cfp.getParams().getDouble("MPLPConvergenceThreshold");
						hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, numMPLPIters, convergenceThreshold);
						order = new StaticScoreHMeanAStarOrder();
					}
					
					// init the A* tree
					ConfAStarTree tree = new ConfAStarTree(order, gscorer, hscorer, rcs);
					tree.initProgress();
					return tree;
				}
			}; 
		}
	}
}
