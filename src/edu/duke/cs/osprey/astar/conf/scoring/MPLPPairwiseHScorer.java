package edu.duke.cs.osprey.astar.conf.scoring;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MPLPPairwiseHScorer implements AStarScorer {
	
	private EnergyMatrix emat;
	
	public MPLPPairwiseHScorer(EnergyMatrix emat) {
		this.emat = emat;
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// run the traditional A* heuristic to get the initial conformation
		
		// init the undefined positions
		int[] undefinedRCs = new int[confIndex.getNumUndefined()];
		Arrays.fill(undefinedRCs, -1);
		
		double energy = 0;
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double minPos1Energy = Double.POSITIVE_INFINITY;
			int minPos1RC = -1;
			
			for (int rc1 : rcs.get(pos1)) {
			
				// undefined single energy
				double pos1Energy = emat.getOneBody(pos1, rc1);
				
				// undefined-defined pairwise energy
				for (int j=0; j<confIndex.getNumDefined(); j++) {
					int pos2 = confIndex.getDefinedPos()[j];
					int rc2 = confIndex.getDefinedRCs()[j];
					if (pos2 < pos1) {
						pos1Energy += emat.getPairwise(pos1, rc1, pos2, rc2);
					}
				}
				
				// undefined-undefined pairwise energy
				for (int j=0; j<confIndex.getNumUndefined(); j++) {
					int pos2 = confIndex.getUndefinedPos()[j];
					
					if (pos2 >= pos1) {
						continue;
					}
					
					double minPos2Energy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						minPos2Energy = Math.min(minPos2Energy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					pos1Energy += minPos2Energy;
				}
				
				if (pos1Energy < minPos1Energy) {
					minPos1Energy = pos1Energy;
					minPos1RC = rc1;
				}
			}
			
			// save the conf
			undefinedRCs[i] = minPos1RC;
			
			energy += minPos1Energy;
		}
		
		// compute the MPLP node-based messages
		// TODO
		
		return energy;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		throw new Error("not implemented");
	}
}
