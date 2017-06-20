package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class NodeUpdater implements MPLPUpdater {
	
	@Override
	public void update(MessageVars lambdas, EnergyMatrix emat) {
		
		// corrected NMPLP (in the notation of the correction paper)
		
		// delta_ji(xi) = -1/(1 + Ni)*gamma_i(xi) + gamma_ji(xi)
		// delta_ij(xj) = -0.5*delta_{j-i}(xj) + 0.5*max_xi [ theta_ij(xi,xj) + 2/(1+Ni)*gamma_i(xi) - gamma_ji(xi) ]
		// where:
		// delta vars are the updatable vars, updated simultaneously for each i
		// delta_{j-i}(xj) = sum_{k in N(j)\i} delta_kj(xj)
		// gamma_ji(xi) = max_xj [ theta_ij(xi,xj) + delta_{j-i}(xj) ]
		// gamma_i(xi) = sum_{j in N(i)} gamma_ji(xi)
		
		// translated into "code notation":
		
		// delta_{posi2,posi1}(rci1) = gamma_{posi2,posi1}(rci1) - gamma_posi1(rci1)/numUndefined
		// delta_{posi1,posi2}(rci2) = [
		//                             max_rci1 [ theta_{posi1,posi2}(rci1,rci2) + 2/numUndefined*gamma_posi1(rci1) - gamma_{posi2,posi1}(rci1) ]
		//                             - delta_{posi2-posi1}(rci2)
		//                           ]/2
		// gamma_{posi2,posi1}(rci1) = max_rci2 [ theta_{posi1,posi2}(rci1,rci2) + delta_{posi2-posi1}(rci2) ]
		// gamma_pos1(rci1) = sum_{posi2 != posi1} gamma_{posi2,posi1}(rci1)
		
		// time complexity
		// g2 = O(n*m)
		// g1 = O(n*g2)
		//    = O(n^2*m)
		// O(n*n*( m*(g2+g1) + m*(m*(g2+g1) + n) ))
		// = O(n^2*( n*m^2 + n^2*m^2 + n*m^3 + n^2*m^3 + m*n ))
		// = O(n^3*m^2 + n^4*m^2 + n^3*m^3 + n^4*m^3 + n^3*m)
		// = O(n^4*m^3)
		
		// with cached gammas:
		// O(n*n*m*n*m) + O(n*n*m*n) + O(n*n*( m + m*m + m*n ))
		// = O(n^3*m^2) + O(n^3*m) + O(n^2*m + n^2*m^2 + n^3*m)
		// = O(n^3*m^2)
		
		// with cached sums (but not cached gammas):
		// g2 = O(m)
		// g1 = O(nm)
		// O(n*n*( m*(g2+g1) + m*m*(g2+g1) ))
		// = O(n^2*( m^2 + n*m^2 + m^3 + n*m^3 ))
		// = O(n^2*m^2 + n^3*m^2 + n^2*m^3 + n^3*m^3)
		// = O(n^3*m^3)
		
		// with cached sums and cached gammas
		// O(n*n*( m + m*m ))
		// = O(n^2*m + n^2*m^2)
		// = O(n^2*m^2)
		// just as fast per update as EMPLP! =)
		
		// NOTE: don't copy-on-write to the lambda vars!
		// the algorithm is designed for immediate updates to the vars
		// as long as the outer loops is over pos1
		
		RCs rcs = lambdas.getRCs();
		ConfIndex confIndex = lambdas.getConfIndex();
		
		for (int posi1=0; posi1<confIndex.numUndefined; posi1++) {
			int pos1 = confIndex.undefinedPos[posi1];
			
			// calculate the gammas
			// NOTE: as far as I know, this precalculation can't be moved outside of the loop =(
			MessageVars gammas = new MessageVars(rcs, confIndex);
			for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
				int pos2 = confIndex.undefinedPos[posi2];
				
				if (pos2 == pos1) {
					continue;
				}
				
				for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
					int rc1 = rcs.get(pos1, rci1);
					
					double minVal = Double.POSITIVE_INFINITY;
					for (int rci2=0; rci2<rcs.getNum(pos2); rci2++) {
						int rc2 = rcs.get(pos2, rci2);
						
						double theta = emat.getPairwise(pos1, rc1, pos2, rc2);
						double delta = lambdas.getEnergyWithout(posi2, rci2, posi1);
						
						minVal = Math.min(minVal, theta + delta);
					}
					
					gammas.set(posi2, posi1, rci1, minVal);
				}
			}
		
			// update the lambdas
			for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
				int pos2 = confIndex.undefinedPos[posi2];
				
				if (pos2 == pos1) {
					continue;
				}
				
				// delta_{posi2,posi1}(rci1) = gamma_{posi2,posi1}(rci1) - gamma_posi1(rci1)/numUndefined
				for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
					
					double gamma = gammas.get(posi2, posi1, rci1);
					double gammaEnergy = gammas.getEnergy(posi1, rci1)/confIndex.numUndefined;
					
					double lambda;
					if (Double.isFinite(gammaEnergy)) {
						lambda = gamma - gammaEnergy;
					} else {
						lambda = Double.POSITIVE_INFINITY;
					}
					
					lambdas.set(posi2, posi1, rci1, lambda);
				}
				
				// delta_{posi1,posi2}(rci2) = [
				//    max_rci1 [
				//       theta_{posi1,posi2}(rci1,rci2)
				//       + 2/numUndefined*gamma_posi1(rci1)
				//       - gamma_{posi2,posi1}(rci1)
				//    ]
				//    - delta_{posi2-posi1}(rci2)
				// ]/2
				for (int rci2=0; rci2<rcs.getNum(pos2); rci2++) {
					int rc2 = rcs.get(pos2, rci2);
					
					double minVal = Double.POSITIVE_INFINITY;
					for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
						int rc1 = rcs.get(pos1, rci1);
						double theta = emat.getPairwise(pos1, rc1, pos2, rc2);
						double gamma1 = gammas.getEnergy(posi1, rci1);
						double gamma2 = gammas.get(posi2, posi1, rci1);
						if (Double.isFinite(theta) && Double.isFinite(gamma1) && Double.isFinite(gamma2)) {
							minVal = Math.min(minVal, theta + 2*gamma1/confIndex.numUndefined - gamma2);
						}
					}
					
					double energyWithout = lambdas.getEnergyWithout(posi2, rci2, posi1);
					if (Double.isFinite(energyWithout)) {
						minVal = (minVal - energyWithout)/2;
					} else {
						minVal = Double.POSITIVE_INFINITY;
					}
					
					lambdas.set(posi1, posi2, rci2, minVal);
				}
			}
		}
	}
}
