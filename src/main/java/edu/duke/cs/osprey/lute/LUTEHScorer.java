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

package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;


public class LUTEHScorer implements AStarScorer {

	public final LUTEConfEnergyCalculator ecalc;

	public LUTEHScorer(LUTEConfEnergyCalculator ecalc) {
		this.ecalc = ecalc;
	}

	@Override
	public AStarScorer make() {
		return new LUTEHScorer(ecalc);
	}


	@Override
	public double calc(ConfIndex index, RCs rcs) {

		// for the moment, LUTE only uses up to triple tuples, so we only need to check up to triples here too

		// TODO: this implementation is very naive... optimize it?
		// although, design algos without minimization are stupidly fast... we probably don't need to optimize this yet?

		double hscore = 0;

		// get the score for each undefined position
		for (int i=0; i<index.numUndefined; i++) {
			int pos1 = index.undefinedPos[i];

			// min over possible assignments to pos1
			double pos1Energy = Double.POSITIVE_INFINITY;
			for (int rc1 : rcs.get(pos1)) {

				double rc1Energy = ecalc.getEnergy(pos1, rc1);

				// interactions with defined residues
				for (int j=0; j<index.numDefined; j++) {
					int pos2 = index.definedPos[j];
					int rc2 = index.definedRCs[j];

					rc1Energy += ecalc.getEnergy(pos1, rc1, pos2, rc2);

					for (int k=0; k<j; k++) {
						int pos3 = index.definedPos[k];
						int rc3 = index.definedRCs[k];

						rc1Energy += ecalc.getEnergy(pos1, rc1, pos2, rc2, pos3, rc3);
					}
				}

				// interactions with undefined residues
				for (int j=0; j<i; j++) {
					int pos2 = index.undefinedPos[j];

					// min over possible assignments to pos2
					double minrc2Energy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {

						// pair with pos2
						double rc2Energy = ecalc.getEnergy(pos1, rc1, pos2, rc2);

						// triples with defined positions
						for (int k=0; k<index.numDefined; k++) {
							int pos3 = index.definedPos[k];
							int rc3 = index.definedRCs[k];
							rc2Energy += ecalc.getEnergy(pos1, rc1, pos2, rc2, pos3, rc3);
						}

						// triples with undefined positions
						for (int k=0; k<j; k++) {
							int pos3 = index.undefinedPos[k];

							// min over rcs
							double minrc3Energy = Double.POSITIVE_INFINITY;
							for (int rc3 : rcs.get(pos3)) {
								double rc3Energy = ecalc.getEnergy(pos1, rc1, pos2, rc2, pos3, rc3);
								minrc3Energy = Math.min(minrc3Energy, rc3Energy);
							}

							rc2Energy += minrc3Energy;
						}

						minrc2Energy = Math.min(minrc2Energy, rc2Energy);
					}

					rc1Energy += minrc2Energy;
				}

				pos1Energy = Math.min(pos1Energy, rc1Energy);
			}

			hscore += pos1Energy;
		}

		return hscore;
	}
}
