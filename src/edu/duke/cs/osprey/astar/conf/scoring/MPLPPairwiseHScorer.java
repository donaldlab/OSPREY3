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
