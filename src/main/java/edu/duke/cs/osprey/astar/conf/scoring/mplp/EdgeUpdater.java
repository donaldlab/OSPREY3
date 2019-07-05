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

package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class EdgeUpdater implements MPLPUpdater {
	
	@Override
	public void update(MessageVars lambdas, EnergyMatrix emat) {
		
		// lambda_ji(xi) = -0.5*lambda_{i-j}(xi) + 0.5*max_xj [ lamda_{j-i}(xj) + theta_ij(xi,xj) ]
		// and with i,j reversed
		
		// lambda_{posi2,posi1}(rci1) = [
		//                              max_{rci2} [ lambda_{posi2-posi1}(rci2) + theta_{posi1,posi2)(rci1,rci2) ]
		//                              - lamba_{posi1-posi2}(rci1)
		//                            ]/2
		// and with posi1,posi2 reversed
		
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
		
		ConfIndex confIndex = lambdas.getConfIndex();
		
		for (int posi1=0; posi1<confIndex.numUndefined; posi1++) {
			int pos1 = confIndex.undefinedPos[posi1];
			
			for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
				int pos2 = confIndex.undefinedPos[posi2];
				
				if (pos2 >= pos1) {
					continue;
				}
				
				update(lambdas, emat, posi1, posi2);
				update(lambdas, emat, posi2, posi1);
			}
		}
	}
	
	private void update(MessageVars lambdas, EnergyMatrix emat, int posi1, int posi2) {

		RCs rcs = lambdas.getRCs();
		ConfIndex confIndex = lambdas.getConfIndex();
		
		int pos1 = confIndex.undefinedPos[posi1];
		int pos2 = confIndex.undefinedPos[posi2];

		for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
			int rc1 = rcs.get(pos1, rci1);
			
			double minEnergy = Double.POSITIVE_INFINITY;
			for (int rci2=0; rci2<rcs.getNum(pos2); rci2++) {
				int rc2 = rcs.get(pos2, rci2);
				double energy = lambdas.getEnergyWithout(posi2, rci2, posi1)
					+ emat.getPairwise(pos1, rc1, pos2, rc2);
				minEnergy = Math.min(minEnergy, energy);
			}
			
			double energyWithout = lambdas.getEnergyWithout(posi1, rci1, posi2);
			if (Double.isFinite(energyWithout)) {
				minEnergy = (minEnergy - energyWithout)/2;
			} else {
				minEnergy = Double.POSITIVE_INFINITY;
			}
			
			lambdas.set(posi2, posi1, rci1, minEnergy);
		}
	}
}
