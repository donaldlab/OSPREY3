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
		
		// delta_{pos2,pos1}(rci1) = gamma_{pos2,pos1}(rci1) - gamma_pos1(rci1)/numUndefined
		// delta_{pos1,pos2}(rci2) = [
		//                             max_rci1 [ theta_{pos1,pos2}(rci1,rci2) + 2/numUndefined*gamma_pos1(rci1) - gamma_{pos2,pos1}(rci1) ]
		//                             - delta_{pos2-pos1}(rci2)
		//                           ]/2
		// gamma_{pos2,pos1}(rci1) = max_rci2 [ theta_{pos1,pos2}(rci1,rci2) + delta_{pos2-pos1}(rci2) ]
		// gamma_pos1(rci1) = sum_{pos2 != pos1} gamma_{pos2,pos1}(rci1)
		
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
		
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			for (int j=0; j<confIndex.getNumUndefined(); j++) {
				int pos2 = confIndex.getUndefinedPos()[j];
				
				if (pos2 == pos1) {
					continue;
				}
				
				// delta_{pos2,pos1}(rci1) = gamma_{pos2,pos1}(rci1) - gamma_pos1(rci1)/numUndefined
				for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
					double delta = nmplpCalcGamma(rcs, confIndex, pos2, pos1, rci1, lambdas, emat)
						- nmplpCalcGamma(rcs, confIndex, pos1, rci1, lambdas, emat)/confIndex.getNumUndefined();
					lambdas.set(pos2, pos1, rci1, delta);
				}
				
				// delta_{pos1,pos2}(rci2) = [
				//    max_rci1 [
				//       theta_{pos1,pos2}(rci1,rci2)
				//       + 2/numUndefined*gamma_pos1(rci1)
				//       - gamma_{pos2,pos1}(rci1)
				//    ]
				//    - delta_{pos2-pos1}(rci2)
				// ]/2
				for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
					int rc2 = rcs.get(pos2)[rci2];
					
					double minVal = Double.POSITIVE_INFINITY;
					for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
						int rc1 = rcs.get(pos1)[rci1];
						double theta = emat.getPairwise(pos1, rc1, pos2, rc2);
						double gamma1 = nmplpCalcGamma(rcs, confIndex, pos1, rci1, lambdas, emat);
						double gamma2 = nmplpCalcGamma(rcs, confIndex, pos2, pos1, rci1, lambdas, emat);
						minVal = Math.min(minVal, theta + 2*gamma1/confIndex.getNumUndefined() - gamma2);
					}
					minVal -= lambdas.getEnergyWithout(pos2, rci2, pos1);
					minVal /= 2;
					lambdas.set(pos1, pos2, rci2, minVal);
				}
			}
		}
	}
	
	private double nmplpCalcGamma(RCs rcs, ConfIndex confIndex, int pos2, int pos1, int rci1, MessageVars lambdas, EnergyMatrix emat) {
		// gamma_{pos2,pos1}(rci1) = max_rci2 [ theta_{pos1,pos2}(rci1,rci2) + delta_{pos2-pos1}(rci2) ]
		int rc1 = rcs.get(pos1)[rci1];
		double minVal = Double.POSITIVE_INFINITY;
		for (int rci2=0; rci2<rcs.get(pos2).length; rci2++) {
			int rc2 = rcs.get(pos2)[rci2];
			double theta = emat.getPairwise(pos1, rc1, pos2, rc2);
			double delta = lambdas.getEnergyWithout(pos2, rci2, pos1);
			minVal = Math.min(minVal, theta + delta);
		}
		return minVal;
	}
	
	private double nmplpCalcGamma(RCs rcs, ConfIndex confIndex, int pos1, int rci1, MessageVars lambdas, EnergyMatrix emat) {
		// gamma_pos1(rci1) = sum_{pos2 != pos1} gamma_{pos2,pos1}(rci1)
		
		// start with single energy
		int rc1 = rcs.get(pos1)[rci1];
		double sum = emat.getOneBody(pos1, rc1);
		
		// add undefined-defined energy
		for (int j=0; j<confIndex.getNumDefined(); j++) {
			int pos2 = confIndex.getDefinedPos()[j];
			int rc2 = confIndex.getDefinedRCs()[j];
			sum += emat.getPairwise(pos1, rc1, pos2, rc2);
		}
		
		// add undefined-undefined energy
		for (int j=0; j<confIndex.getNumUndefined(); j++) {
			int pos2 = confIndex.getUndefinedPos()[j];
			if (pos2 != pos1) {
				sum += nmplpCalcGamma(rcs, confIndex, pos2, pos1, rci1, lambdas, emat);
			}
		}
		
		return sum;
	}
}
