package edu.duke.cs.osprey.astar.conf.ranking;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;


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

	private GScoreAStarNode lowerNode = new GScoreAStarNode();
	private GScoreAStarNode upperNode = new GScoreAStarNode();
	private int[] oneRc = { 0 };

	@Override
	public int getNextPosition(ConfRanker ranker, ConfIndex confIndex, RCs rcs, double queryScore) {

		// get g-scores for this sub-tree
		lowerNode.gscore = ranker.lowerGScorer.calc(confIndex, rcs);
		upperNode.gscore = ranker.upperGScorer.calc(confIndex, rcs);

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

				// backup the conf index's node
				ConfAStarNode node = confIndex.node;

				// approximate the optimal sub-tree min,max scores using the A* heuristic
				// then see if any of the sub-sub-trees could be pruned
				confIndex.node = lowerNode;
				double minScore = 0
					+ ranker.lowerGScorer.calcDifferential(confIndex, subRCs, pos, rc)
					+ ranker.lowerHScorer.calcDifferential(confIndex, subRCs, pos, rc);
				if (minScore > queryScore) {
					posScore += fullRCScore;
				} else {

					confIndex.node = upperNode;
					double maxScore = 0
						- ranker.upperGScorer.calcDifferential(confIndex, subRCs, pos, rc)
						- ranker.upperHScorer.calcDifferential(confIndex, subRCs, pos, rc);
					if (maxScore <= queryScore) {
						posScore += fullRCScore;
					}
				}

				// restore the conf index state
				confIndex.node = node;
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
