package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class EdgeUpdater implements MPLPUpdater {
	
	@Override
	public void update(MessageVars lambdas, EnergyMatrix emat) {
		
		// lambda_ji(xi) = -0.5*lambda_{i-j}(xi) + 0.5*max_xj [ lamda_{j-i}(xj) + theta_ij(xi,xj) ]
		// and with i,j reversed
		
		// lambda_{pos2,pos1}(rci1) = [
		//                              max_{rci2} [ lambda_{pos2-pos1}(rci2) + theta_{pos1,pos2)(rci1,rci2) ]
		//                              - lamba_{pos1-pos2}(rci1)
		//                            ]/2
		// and with pos1,pos2 reversed
		
		// time complexity
		// O(n*n*2*m*(m*n + n))
		//  = O(2n^2*m*(nm + n))
		//  = O(2n^3*m^2 + 2n^3*m)
		
		// with caching of lambda sums
		// O(n*n*2*m( m ))
		//  = O(2*n^2*m^2)
		
		// NOTE: don't copy-on-write to the lambda vars!
		// the algorithm is designed for immediate updates to the lambda vars
		// as long as the outer loops are over pos1,pos2
		
		RCs rcs = lambdas.getRCs();
		ConfIndex confIndex = lambdas.getConfIndex();
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos2 >= pos1) {
					continue;
				}
				
				for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
					int rc1 = rcs.get(pos1)[rci1];
					
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
						int rc2 = rcs.get(pos2)[rci2];
						double energy = lambdas.getEnergyWithout(pos2, rci2, pos1)
							+ emat.getPairwise(pos1, rc1, pos2, rc2);
						minEnergy = Math.min(minEnergy, energy);
					}
					
					minEnergy = (minEnergy - lambdas.getEnergyWithout(pos1, rci1, pos2))/2;
					
					lambdas.set(pos2, pos1, rci1, minEnergy);
				}
				
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					int rc2 = rcs.get(pos2)[rci2];
					
					double minEnergy = Double.POSITIVE_INFINITY;
					for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
						int rc1 = rcs.get(pos1)[rci1];
						double energy = lambdas.getEnergyWithout(pos1, rci1, pos2)
							+ emat.getPairwise(pos1, rc1, pos2, rc2);
						minEnergy = Math.min(minEnergy, energy);
					}
					
					minEnergy = (minEnergy - lambdas.getEnergyWithout(pos2, rci2, pos1))/2;
					
					lambdas.set(pos1, pos2, rci2, minEnergy);
				}
			}
		}
	}
}
