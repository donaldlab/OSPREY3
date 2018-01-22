package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.Comparator;


public class StaticOptimalPruningOrderer implements ConfRanker.Orderer {

	private int[] order = null;

	@Override
	public int getNextPosition(ConfRanker ranker, ConfIndex confIndex, RCs rcs, double queryScore) {

		// compute the order the first time if needed
		if (order == null) {
			computeOrder(ranker, queryScore);
		}

		// return the first unassigned pos in the order
		for (int pos : order) {
			if (confIndex.isUndefined(pos)) {
				return pos;
			}
		}

		// this shouldn't be possible
		throw new Error("unpossible?");
	}

	private void computeOrder(ConfRanker ranker, double queryScore) {
		order = ranker.confSpace.positions.stream()
			.sorted(Comparator.comparing((SimpleConfSpace.Position pos) -> {

				// score the position by the number of confs pruned

				BigInteger numConfsPruned = BigInteger.ZERO;

				for (int rc : ranker.rcs.get(pos.index)) {

					int [] subConfMask = ranker.noAssignmentsMask.clone();
					subConfMask[pos.index] = rc;
					RCs subRCs = ranker.makeRCs(subConfMask);

					// can the confs in this sub-tree can be pruned?
					if (ranker.getMinScore(subRCs) > queryScore) {
						numConfsPruned = numConfsPruned.add(subRCs.getNumConformations());
					} else if (ranker.getMaxScore(subRCs) <= queryScore) {
						numConfsPruned = numConfsPruned.add(subRCs.getNumConformations());
					}
				}

				return numConfsPruned;
			}).reversed())
			.mapToInt((pos) -> pos.index)
			.toArray();
	}
}