package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.math.BigInteger;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;


public class StaticOptimalPruningOrderer implements ConfRanker.Orderer {

	private List<SimpleConfSpace.Position> order = null;

	@Override
	public SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore) {

		// compute the order the first time if needed
		if (order == null) {
			computeOrder(ranker, queryScore);
		}

		// return the first unassigned pos in the order
		for (SimpleConfSpace.Position pos : order) {
			if (unassignedPositions.contains(pos)) {
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
					RCs subRCs = ranker.makeSubRCs(subConfMask);

					// can the confs in this sub-tree can be pruned?
					double minScore = ranker.getMinScore(subRCs);
					if (minScore > queryScore) {
						numConfsPruned = numConfsPruned.add(subRCs.getNumConformations());
						continue;
					}

					double maxScore = ranker.getMaxScore(subRCs);
					if (maxScore <= queryScore) {
						numConfsPruned = numConfsPruned.add(subRCs.getNumConformations());
					}
				}

				return numConfsPruned;
			}).reversed())
			.collect(Collectors.toList());
	}
}