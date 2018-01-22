package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.List;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;

/**
 * chooses the next position that optimizes conf sub-tree pruning,
 * but is expensive to evaluate in practice
 */
public class DynamicHeuristicPruningOrderer implements ConfRanker.Orderer {

	class GScoreAStarNode implements ConfAStarNode {

		public double gscore;

		@Override
		public ConfAStarNode assign(int pos, int rc) {
			throw new UnsupportedOperationException();
		}

		@Override
		public double getGScore() {
			return gscore;
		}

		@Override
		public void setGScore(double val) {
			throw new UnsupportedOperationException();
		}

		@Override
		public double getHScore() {
			throw new UnsupportedOperationException();
		}

		@Override
		public void setHScore(double val) {
			throw new UnsupportedOperationException();
		}

		@Override
		public int getLevel() {
			throw new UnsupportedOperationException();
		}

		@Override
		public void getConf(int[] conf) {
			throw new UnsupportedOperationException();
		}

		@Override
		public void index(ConfIndex index) {
			throw new UnsupportedOperationException();
		}
	}

	private ConfIndex confIndex = null;
	private GScoreAStarNode lowerNode = new GScoreAStarNode();
	private GScoreAStarNode upperNode = new GScoreAStarNode();
	private AStarScorer lowerGScorer = null;
	private AStarScorer lowerHScorer = null;
	private AStarScorer upperGScorer = null;
	private AStarScorer upperHScorer = null;
	private int[] oneRc = { 0 };

	@Override
	public SimpleConfSpace.Position getNextPosition(ConfRanker ranker, int[] confMask, List<SimpleConfSpace.Position> unassignedPositions, double queryScore) {

		// first time init, if needed
		if (confIndex == null) {
			confIndex = new ConfIndex(ranker.confSpace.positions.size());

			RCs rcs = new RCs(ranker.confSpace);

			ConfAStarTree lowerAStar = ranker.makeAStar(ranker.emat, rcs);
			lowerGScorer = lowerAStar.gscorer;
			lowerHScorer = lowerAStar.hscorer;

			ConfAStarTree upperAStar = ranker.makeAStar(ranker.negatedEmat, rcs);
			upperGScorer = upperAStar.gscorer;
			upperHScorer = upperAStar.hscorer;
		}

		// populate the conf index with the mask
		confIndex.numDefined = 0;
		confIndex.numUndefined = 0;
		for (SimpleConfSpace.Position pos : ranker.confSpace.positions) {
			int rc = confMask[pos.index];
			if (rc == -1) {
				confIndex.undefinedPos[confIndex.numUndefined] = pos.index;
				confIndex.numUndefined++;
			} else {
				confIndex.definedPos[confIndex.numDefined] = pos.index;
				confIndex.definedRCs[confIndex.numDefined] = confMask[pos.index];
				confIndex.numDefined++;
			}
		}

		// get g-scores for this sub-tree
		RCs rcs = ranker.makeSubRCs(confMask);
		lowerNode.gscore = lowerGScorer.calc(confIndex, rcs);
		upperNode.gscore = upperGScorer.calc(confIndex, rcs);

		// allocate space for the sub-rcs
		RCs subRCs = new RCs(rcs);

		SimpleConfSpace.Position nextPos = null;
		double bestPosScore = Double.NEGATIVE_INFINITY;
		for (SimpleConfSpace.Position unassignedPos : unassignedPositions) {

			double posScore = 0.0;
			int[] posRCs = rcs.get(unassignedPos.index);

			for (int rc : posRCs) {

				// update the sub-RCs with this pos and rc
				for (int pos=0; pos<rcs.getNumPos(); pos++) {
					if (pos == unassignedPos.index) {
						subRCs.set(pos, oneRc);
					} else {
						subRCs.set(pos, rcs.get(pos));
					}
				}
				oneRc[0] = rc;

				double fullRCScore = 1.0/posRCs.length;

				// approximate the optimal sub-tree min,max scores using the A* heuristic
				// then see if any of the sub-sub-trees could be pruned
				confIndex.node = lowerNode;
				double minScore = 0
					+ lowerGScorer.calcDifferential(confIndex, subRCs, unassignedPos.index, rc)
					+ lowerHScorer.calcDifferential(confIndex, subRCs, unassignedPos.index, rc);
				if (minScore > queryScore) {
					posScore += fullRCScore;
					continue;
				}

				confIndex.node = upperNode;
				double maxScore = 0
					- upperGScorer.calcDifferential(confIndex, subRCs, unassignedPos.index, rc)
					- upperHScorer.calcDifferential(confIndex, subRCs, unassignedPos.index, rc);
				if (maxScore <= queryScore) {
					posScore += fullRCScore;
					continue;
				}

				/* NOPE, this doesn't work well
				// can't prune anything, try to award a partial score for closeness of the query score to the bounds
				double centerScore = (minScore + maxScore)/2.0;
				double halfWidth = (maxScore - minScore)/2.0;
				double distToCenter = Math.abs(queryScore - centerScore);
				double distRatio = distToCenter/halfWidth;
				posScore += distRatio*fullRCScore;
				*/

				//log("\tconfMask: %s  next %d,%d   bounds [%.4f,%.4f]   query: %.4f   dist: %.4f   ratio: %.12f", Arrays.toString(confMask), unassignedPos.index, rc, minScore, maxScore, queryScore, distToCenter, distRatio);
			}

			// update the best position based on ratio of sub-trees pruned
			if (posScore > bestPosScore) {
				bestPosScore = posScore;
				nextPos = unassignedPos;
			}

			//log("pos %d   sub-trees pruned: %d/%d   pos score: %.4f", unassignedPos.index, numSubTreesPruned, rcs.get(unassignedPos.index).length, posScore);
		}
		assert (nextPos != null);

		if (bestPosScore == 0.0) {

			// no pruning possible, order to minimize number of sub-trees instead
			for (SimpleConfSpace.Position unassignedPos : unassignedPositions) {

				double posScore = -rcs.getNum(unassignedPos.index);

				if (posScore > bestPosScore) {
					bestPosScore = posScore;
					nextPos = unassignedPos;
				}
			}
		}

		//log("best post to try next: %d -> %s",nextPos.index, formatBig(bestNextPosScore));

		return nextPos;
	}
}
