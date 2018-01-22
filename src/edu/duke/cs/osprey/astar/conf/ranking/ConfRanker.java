package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.math.BigInteger;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;

public class ConfRanker {

	public static interface ConfAStarConfigurer {
		void configure(ConfAStarTree.Builder builder);
	}

	public static interface Orderer {
		SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore);
	}

	private class BigProgress {

		public final BigInteger total;
		public BigInteger below = BigInteger.ZERO;
		public BigInteger above = BigInteger.ZERO;

		public final long startTimeNs = System.nanoTime();
		public long reportIntervalMs = 5000;

		private long lastReportTimeMs = 0;

		public BigProgress(BigInteger total) {
			this.total = total;
		}

		public void incrementBelow() {
			incrementBelow(BigInteger.ONE);
		}

		public void incrementAbove() {
			incrementAbove(BigInteger.ONE);
		}

		public void incrementBelow(BigInteger val) {
			below = below.add(val);
			writeReportIfNeeded();
		}

		public void incrementAbove(BigInteger val) {
			above = above.add(val);
			writeReportIfNeeded();
		}

		public void writeReportIfNeeded() {

			// should we write a report?
			if (reportProgress) {
				long timeMs = System.currentTimeMillis();
				if (timeMs - lastReportTimeMs >= reportIntervalMs) {
					lastReportTimeMs = timeMs;
					System.out.println(getReport());
				}
			}
		}

		public String getReport() {
			return String.format("progress: [%e,%e] of %e  (%.6f%%)   %s",
				below.doubleValue(),
				total.doubleValue() - above.doubleValue(),
				total.doubleValue(),
				below.add(above).doubleValue()/total.doubleValue()*100.0,
				TimeFormatter.format(System.nanoTime() - startTimeNs, 2)
			);
		}
	}

	private class Entry {

		public final int[] confMask;
		public final int numAssignments;

		public Entry(int[] confMask) {

			this.confMask = confMask;

			int count = 0;
			for (int rc : confMask) {
				if (rc != NotAssigned) {
					count++;
				}
			}
			this.numAssignments = count;
		}

		public Entry(Entry parent, SimpleConfSpace.Position pos, int rc) {
			this(assign(parent.confMask, pos, rc));
		}
	}

	private static int[] assign(int[] confMask, SimpleConfSpace.Position pos, int rc) {
		int[] subConfMask = confMask.clone();
		subConfMask[pos.index] = rc;
		return subConfMask;
	}

	private static final int NotAssigned = -1;

	public final SimpleConfSpace confSpace;
	public final RCs rcs;
	public final EnergyMatrix emat;
	public final Orderer orderer;
	public final ConfAStarConfigurer astarConfigurer;

	public final int[] noAssignmentsMask;

	public boolean reportProgress = false;

	public final EnergyMatrix negatedEmat;
	private final AStarScorer scorer;
	private final ConfIndex confIndex;
	private final PriorityQueue<Entry> queue;

	public ConfRanker(SimpleConfSpace confSpace, RCs rcs, EnergyMatrix emat, Orderer orderer, ConfAStarConfigurer astarConfigurer) {

		this.confSpace = confSpace;
		this.rcs = rcs;
		this.emat = emat;
		this.orderer = orderer;
		this.astarConfigurer = astarConfigurer;

		negatedEmat = new NegatedEnergyMatrix(confSpace, emat);

		// make the all-unassigned conf mask
		noAssignmentsMask = new int[confSpace.positions.size()];
		Arrays.fill(noAssignmentsMask, NotAssigned);

		// get the A* scorer
		scorer = makeAStar(emat, rcs).gscorer;

		// allocate the conf index for conf scoring, and prep for fully-assigned confs
		confIndex = new ConfIndex(confSpace.positions.size());
		confIndex.numDefined = confSpace.positions.size();
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			 confIndex.definedPos[pos.index] = pos.index;
		}
		confIndex.numUndefined = 0;

		queue = new PriorityQueue<>(Comparator.comparing((entry) -> entry.numAssignments));
	}

	public RCs makeSubRCs(int[] confMask) {
		return new RCs(rcs, (pos, rc) -> {
			return confMask[pos] == -1 || confMask[pos] == rc;
		});
	}

	public ConfAStarTree makeAStar(EnergyMatrix emat, RCs rcs) {
		ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs);
		astarConfigurer.configure(builder);
		return builder.build();
	}

	public double getMinScore(RCs rcs) {
		return makeAStar(emat, rcs)
			.nextConf()
			.getScore();
	}

	public double getMaxScore(RCs rcs) {
		return -makeAStar(negatedEmat, rcs)
			.nextConf()
			.getScore();
	}

	public BigInteger getNumConfsAtMost(double queryScore) {

		BigProgress progress = new BigProgress(rcs.getNumConformations());

		// NOTE: conf masks are always sequences of assigned positions, then unassigned positions
		// e.g., -1, -1, -1
		// or 0, 5, -1
		// but never -1, -1, 2

		// NOTE: use a priority traversal of the tree, so we can stop early and still maybe get useful bounds on the rank
		// DFS spends a lot of time up front making only minor adjustments to the counts and is not good for early stopping

		// start with all unassigned
		queue.add(new Entry(noAssignmentsMask));

		while (!queue.isEmpty()) {

			Entry entry = queue.poll();

			// get the unassigned positions, if any
			List<SimpleConfSpace.Position> unassignedPositions = confSpace.positions.stream()
				.filter((pos) -> entry.confMask[pos.index] == -1)
				.collect(Collectors.toList());

			// if all positions are assigned (and confMask is really just a conf),
			// then just compare conf scores directly
			if (unassignedPositions.isEmpty()) {

				// update the conf index with the conf RCs
				for (SimpleConfSpace.Position pos : confSpace.positions) {
					confIndex.definedRCs[pos.index] = entry.confMask[pos.index];
				}
				double confScore = scorer.calc(confIndex, rcs);

				if (confScore <= queryScore) {
					progress.incrementBelow();
				} else {
					progress.incrementAbove();
				}

				// leaf node, we're done here
				continue;
			}

			// get the RCs for this sub-tree using the conf mask
			RCs rcs = makeSubRCs(entry.confMask);

			// no flexibility? no conformations
			// (this probably shouldn't happen, but I guess it can't hurt if it does)
			if (!rcs.hasConfs()) {
				continue;
			}

			// is the whole sub-tree above the query?
			if (getMinScore(rcs) > queryScore) {

				progress.incrementAbove(rcs.getNumConformations());

				// sub-tree pruned
				continue;
			}

			// is the whole sub-tree below the query?
			if (getMaxScore(rcs) <= queryScore) {

				progress.incrementBelow(rcs.getNumConformations());

				// sub-tree pruned
				continue;
			}

			// didn't prune anything yet, but update progress if needed
			progress.writeReportIfNeeded();

			// what's the next position we should assign?
			SimpleConfSpace.Position nextPos = orderer.getNextPosition(this, entry.confMask, unassignedPositions, queryScore);

			// the sub-tree is part below and part above, so queue up the sub-sub-trees
			for (int rc : rcs.get(nextPos.index)) {
				queue.add(new Entry(entry, nextPos, rc));
			}
		}

		return progress.below;
	}
}
