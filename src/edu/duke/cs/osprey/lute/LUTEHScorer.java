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
