package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.List;

/**
 * chooses the next position that optimizes conf sub-tree pruning,
 * but is expensive to evaluate in practice
 */
public class DynamicOptimalPruningOrderer implements ConfRanker.Orderer {

	@Override
	public SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore) {

		// TODO: could do in parallel?
		SimpleConfSpace.Position nextPos = null;
		BigInteger bestNextPosScore = null;
		for (SimpleConfSpace.Position unassignedPos : unassignedPositions) {

			int numSubTreesPruned = 0;

			// NOTE: all assignments at this pos have sub-trees with the same number of confs
			BigInteger numConfs = null;

			for (int nextRc : ranker.rcs.get(unassignedPos.index)) {

				int [] subConfMask = confMask.clone();
				subConfMask[unassignedPos.index] = nextRc;
				RCs subRCs = ranker.makeSubRCs(subConfMask);

				if (numConfs == null) {
					numConfs = subRCs.getNumConformations();
				}

				// can the confs in this sub-tree can be pruned?
				double minScore = ranker.getMinScore(subRCs);
				if (minScore > queryScore) {
					numSubTreesPruned++;
					//continue;
				}

				double maxScore = ranker.getMaxScore(subRCs);
				if (maxScore <= queryScore) {
					numSubTreesPruned++;
				}

				//log("\tconfMask: %s   bounds [%.4f,%.4f]  %s", Arrays.toString(subConfMask), minScore, maxScore, isPruned ? "PRUNED" : "");
			}
			assert (numConfs != null);

			// update the best position based on the number of confs pruned
			BigInteger numConfsPruned = numConfs.multiply(BigInteger.valueOf(numSubTreesPruned));
			if (bestNextPosScore == null || numConfsPruned.compareTo(bestNextPosScore) > 0) {
				bestNextPosScore = numConfsPruned;
				nextPos = unassignedPos;
			}

			//log("pos %d   sub-trees pruned: %d   confs pruned: %s", unassignedPos.index, numSubTreesPruned, formatBig(numConfsPruned));
		}
		assert (nextPos != null);

		//log("best post to try next: %d -> %s",nextPos.index, formatBig(bestNextPosScore));

		return nextPos;
	}
}
