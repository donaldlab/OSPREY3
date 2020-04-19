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
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class TriplewiseHScorer implements AStarScorer {

	// TODO: this needs to be heavily optimized before it can ever be useful!

	private EnergyMatrix emat;

	public TriplewiseHScorer(EnergyMatrix emat) {
		this.emat = emat;
	}
	
	@Override
	public TriplewiseHScorer make() {
		return new TriplewiseHScorer(emat);
	}
	
	@Override
	public double calc(ConfIndex index, RCs rcs) {

		double hscore = 0;

		// allocate just one tuple, but update it inside the loops
		RCTuple triple = new RCTuple(0, 0, 0, 0, 0, 0);
		
		// get the score for each undefined position
		for (int i=0; i<index.numUndefined; i++) {
			int pos1 = index.undefinedPos[i];

			// min over possible assignments to pos1
			double minPos1Energy = Double.POSITIVE_INFINITY;
			for (int rc1 : rcs.get(pos1)) {
				
				double rc1Energy = emat.getOneBody(pos1, rc1);

				// interactions with defined residues
				for (int j=0; j<index.numDefined; j++) {
					int pos2 = index.definedPos[j];
					int rc2 = index.definedRCs[j];

					rc1Energy += emat.getPairwise(pos1, rc1, pos2, rc2);

					// add triples, if any
					for (int k=0; k<j; k++) {
						int pos3 = index.definedPos[k];
						if (pos3 >= pos2) {
							break;
						}
						int rc3 = index.definedRCs[k];
						rc1Energy += getTriplewise(emat, triple, pos1, rc1, pos2, rc2, pos3, rc3);
					}
				}

				// interactions with undefined residues
				for (int j=0; j<index.numUndefined; j++) {
					int pos2 = index.undefinedPos[j];
					if (pos2 >= pos1) {
						break;
					}

					// min over possible assignments to pos2
					double minPos2Energy = Double.POSITIVE_INFINITY;
					for (int rc2 : rcs.get(pos2)) {

						double rc2Energy = emat.getPairwise(pos1, rc1, pos2, rc2);

						// add triples, if any
						for (int k=0; k<j; k++) {
							int pos3 = index.definedPos[k];
							if (pos3 >= pos2) {
								break;
							}
							int rc3 = index.definedRCs[k];
							rc2Energy += getTriplewise(emat, triple, pos1, rc1, pos2, rc2, pos3, rc3);
						}

						// TODO: look at a third position to get the full triple-wise effect?

						minPos2Energy = Math.min(minPos2Energy, rc2Energy);
					}

					rc1Energy += minPos2Energy;
				}

				minPos1Energy = Math.min(minPos1Energy, rc1Energy);
			}
		
			hscore += minPos1Energy;
		}
		
		return hscore;
	}

	private double getTriplewise(EnergyMatrix emat, RCTuple triple, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {

		assert (pos1 != pos2);
		assert (pos1 != pos3);
		assert (pos2 != pos3);

		// sort the positions so: pos1 < pos2 < pos3
		// just use a simple swap chain
		if (pos1 > pos2) {
			int swapPos = pos2;
			int swapRc = rc2;
			pos2 = pos1;
			rc2 = rc1;
			pos1 = swapPos;
			rc1 = swapRc;
		}

		if (pos2 > pos3) {
			int swapPos = pos3;
			int swapRc = rc3;
			pos3 = pos2;
			rc3 = rc2;
			pos2 = swapPos;
			rc2 = swapRc;
		}

		if (pos1 > pos2) {
			int swapPos = pos2;
			int swapRc = rc2;
			pos2 = pos1;
			rc2 = rc1;
			pos1 = swapPos;
			rc1 = swapRc;
		}

		// the compiler is apparently smart enough to know this is correct. nice!
		//assert (pos1 < pos2);
		//assert (pos2 < pos3);

		triple.pos.set(0, pos1);
		triple.pos.set(1, pos2);
		triple.pos.set(2, pos3);

		triple.RCs.set(0, rc1);
		triple.RCs.set(1, rc2);
		triple.RCs.set(2, rc3);

		Double energy = emat.getTuple(triple);
		if (energy == null) {
			return 0.0;
		}
		return energy;
	}
}
