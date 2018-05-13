/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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
	public MPLPPairwiseHScorer make() {
		return new MPLPPairwiseHScorer(updater, emat, maxNumIterations, epsilon);
	}

	@Override
	public double calc(ConfIndex confIndex, RCs rcs) {
		
		// init lambdas using the traditional A* heuristic
		// NOTE: we must use these initial values for early stopping to be sound
		MessageVars lambdas = new MessageVars(rcs, confIndex);
		lambdas.initTraditionalAStar(emat);
		
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
}
