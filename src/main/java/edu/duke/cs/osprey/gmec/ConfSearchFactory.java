/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.gmec;

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
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.multistatekstar.MSSearchProblem;
import edu.duke.cs.osprey.multistatekstar.MultiSequenceConfTree;
import edu.duke.cs.osprey.parallelism.Parallelism;
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
					
					if (search.searchNeedsHigherOrderTerms() || search.useEPIC || cfp.hasGMECMutFile()) {
				
						// if we need higher-order or EPIC terms, use the old A* code
						return ConfTree.makeFull(search, pmat, cfp.parseGMECMutFile(search.confSpace));

					} else if (search instanceof MSSearchProblem && !((MSSearchProblem)search).isFullyAssigned()) {
						// we need a multi-sequence conf space
						return new MultiSequenceConfTree((MSSearchProblem)search, emat, pmat);
					}
					
					// when we don't need higher order terms, we can do fast pairwise-only things
					
					// get the appropriate configuration for A*
					AStarScorer gscorer = new PairwiseGScorer(emat);
					AStarScorer hscorer;
					AStarOrder order;
					RCs rcs = new RCs(pmat);
					
					// how many iterations of MPLP should we do?
					int numMPLPIters = cfp.params.getInt("NumMPLPIters");
					if (numMPLPIters <= 0) {
						
						// zero MPLP iterations is exactly the traditional heuristic, so use the fast implementation
						hscorer = new TraditionalPairwiseHScorer(emat, rcs);
						order = new DynamicHMeanAStarOrder();
						
					} else {
						
						// when reference energies are used (and they are by default), edge-based MPLP is waaay
						// faster than node-based, so use edge-based by default too.
						MPLPUpdater updater;
						if (cfp.params.getValue("MPLPAlg").equalsIgnoreCase("node")) {
							updater = new NodeUpdater();
						} else {
							updater = new EdgeUpdater();
						}
						
						double convergenceThreshold = cfp.params.getDouble("MPLPConvergenceThreshold");
						hscorer = new MPLPPairwiseHScorer(updater, emat, numMPLPIters, convergenceThreshold);
						
						// also, always use a static order with MPLP
						// MPLP isn't optimized to do differential node scoring quickly so DynamicHMean is super slow!
						order = new StaticScoreHMeanAStarOrder();
					}
					
					// init the A* tree
					ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, rcs)
						.setCustom(order, gscorer, hscorer)
						.build();
					tree.setParallelism(Parallelism.makeCpu(cfp.params.getInt("AStarThreads")));
					tree.initProgress();
					return tree;
				}
			}; 
		}
	}
}
