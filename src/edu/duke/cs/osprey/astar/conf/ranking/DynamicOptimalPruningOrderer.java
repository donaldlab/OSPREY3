package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;


/**
 * chooses the next position that optimizes conf sub-tree pruning,
 * but is expensive to evaluate in practice
 */
public class DynamicOptimalPruningOrderer implements ConfRanker.Orderer {

	private int[] oneRc = { 0 };

	@Override
	public int getNextPosition(ConfRanker ranker, ConfIndex confIndex, RCs rcs, double queryScore) {

		// allocate space for the sub-rcs
		RCs subRCs = new RCs(rcs);

		int nextPos = -1;
		double bestPosScore = Double.NEGATIVE_INFINITY;

		for (int i=0; i<confIndex.numUndefined; i++) {
			int pos = confIndex.undefinedPos[i];

			double posScore = 0.0;
			int[] posRCs = rcs.get(pos);

			for (int rc : posRCs) {

				// update the sub-RCs with this pos and rc
				for (int pos2=0; pos2<rcs.getNumPos(); pos2++) {
					if (pos2 == pos) {
						subRCs.set(pos2, oneRc);
					} else {
						subRCs.set(pos2, rcs.get(pos2));
					}
				}
				oneRc[0] = rc;

				double fullRCScore = 1.0/posRCs.length;

				// can the confs in this sub-tree be pruned?
				if (ranker.getMinScore(subRCs) > queryScore) {
					posScore += fullRCScore;
				} else if (ranker.getMaxScore(subRCs) <= queryScore) {
					posScore += fullRCScore;
				}
			}

			// update the best position based on ratio of sub-trees pruned
			if (posScore > bestPosScore) {
				bestPosScore = posScore;
				nextPos = pos;
			}
		}
		assert (nextPos >= 0);

		if (bestPosScore == 0.0) {

			// no pruning possible, order to minimize number of sub-trees instead
			for (int i=0; i<confIndex.numUndefined; i++) {
				int pos = confIndex.undefinedPos[i];

				double posScore = -rcs.getNum(pos);

				if (posScore > bestPosScore) {
					bestPosScore = posScore;
					nextPos = pos;
				}
			}
		}

		return nextPos;
	}
}
