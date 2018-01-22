package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import org.jetbrains.annotations.NotNull;

import java.math.BigInteger;
import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class CachedOptimalPruningOrderer implements ConfRanker.Orderer {

	private class PosIndex {

		public static final int Unassigned = -1;

		public int[] values;

		public PosIndex(int[] confMask) {
			values = confMask.clone();
		}

		@Override
		public int hashCode() {
			return Arrays.hashCode(values);
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof PosIndex && Arrays.equals(this.values, ((PosIndex)other).values);
		}

		public boolean isAssigned(SimpleConfSpace.Position pos) {
			return values[pos.index] != Unassigned;
		}

		public void unassign(SimpleConfSpace.Position pos) {
			values[pos.index] = Unassigned;
		}

		@Override
		public String toString() {
			return Arrays.toString(values);
		}

		public int getNumUnassigned() {
			int count = 0;
			for (int val : values) {
				if (val == Unassigned) {
					count++;
				}
			}
			return count;
		}
	}

	private Map<PosIndex,List<SimpleConfSpace.Position>> orders = new HashMap<>();

	@Override
	public SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore) {

		// TEMP
		final int minNumUnassignedPositions = ranker.confSpace.positions.size() - 1;

		PosIndex posIndex = new PosIndex(confMask);

		List<SimpleConfSpace.Position> order = null;
		if (unassignedPositions.size() < minNumUnassignedPositions) {

			// don't want to store too many orders,
			// so if we only have a few unassigned positions left,
			// pick another related order by arbitrarily unassigning positions
			// TODO: find a way to bubble up to the parent order?
			for (SimpleConfSpace.Position pos : ranker.confSpace.positions) {
				if (posIndex.isAssigned(pos)) {
					posIndex.unassign(pos);

					order = orders.get(posIndex);
					if (order != null) {
						break;
					}
				}
			}

		} else {

			// compute the order the first time if needed
			order = orders.get(posIndex);
			if (order == null) {
				// TEMP
				log("computing order for %s", posIndex);
				order = computeOrder(ranker, queryScore, unassignedPositions);
				orders.put(posIndex, order);
			}
		}
		assert (order != null);

		// return the first unassigned pos in the order
		for (SimpleConfSpace.Position pos : order) {
			if (unassignedPositions.contains(pos)) {
				return pos;
			}
		}

		// this shouldn't be possible
		throw new Error("unpossible?");
	}

	private List<SimpleConfSpace.Position> computeOrder(ConfRanker ranker, double queryScore, List<SimpleConfSpace.Position> unassignedPositions) {
		return unassignedPositions.stream()
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