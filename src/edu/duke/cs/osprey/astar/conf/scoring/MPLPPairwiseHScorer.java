package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MessageVars;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MPLPPairwiseHScorer implements AStarScorer {
	
	private MPLPUpdater updater;
	private EnergyMatrix emat;
	private int maxNumIterations;
	private double epsilon;
	
	public MPLPPairwiseHScorer(MPLPUpdater updater, EnergyMatrix emat, int maxNumIterations, double epsilon) {
		this.updater = updater;
		this.emat = emat;
		this.maxNumIterations = maxNumIterations;
		this.epsilon = epsilon;
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// init lambdas using the traditional A* heuristic
		// NOTE: we must use these initial values for early stopping to be sound
		MessageVars lambdas = new MessageVars(rcs, confIndex, emat);
		for (int posi1=0; posi1<confIndex.getNumUndefined(); posi1++) {
			int pos1 = confIndex.getUndefinedPos()[posi1];
			
			for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
				int rc1 = rcs.get(pos1, rci1);
			
				for (int posi2=0; posi2<confIndex.getNumUndefined(); posi2++) {
					int pos2 = confIndex.getUndefinedPos()[posi2];
					
					if (pos2 >= pos1) {
						continue;
					}
					
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {
						minEnergy = Math.min(minEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
					}
					
					lambdas.set(posi2, posi1, rci1, minEnergy);
				}
			}
		}
		
		// run MPLP
		double energy = lambdas.getTotalEnergy();
		for (int i=0; i<maxNumIterations; i++) {
			updater.update(lambdas, emat);
			double newEnergy = lambdas.getTotalEnergy();
			if (Math.abs(newEnergy - energy) < epsilon) {
				break;
			}
			energy = newEnergy;
		}
		return energy;
	}

	@Override
	public double calcDifferential(ConfIndex confIndex, RCs rcs, int nextPos, int nextRc) {
		
		// TODO: implement me
		
		// TEMP: punt to calc()
		return calc(new ConfIndex(confIndex, nextPos, nextRc), rcs);
	}
}
