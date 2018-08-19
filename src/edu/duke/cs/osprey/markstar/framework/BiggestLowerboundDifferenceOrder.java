package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.tools.MathTools;

public class BiggestLowerboundDifferenceOrder implements edu.duke.cs.osprey.astar.conf.order.AStarOrder {

	public final MathTools.Optimizer optimizer;

	public BiggestLowerboundDifferenceOrder() {
		this(MathTools.Optimizer.Maximize);
	}

	public BiggestLowerboundDifferenceOrder(MathTools.Optimizer optimizer) {
		this.optimizer = optimizer;
	}

	private AStarScorer gscorer;
	private AStarScorer hscorer;

	@Override
	public void setScorers(AStarScorer gscorer, AStarScorer hscorer) {
		this.gscorer = gscorer;
		this.hscorer = hscorer;
	}

	@Override
	public boolean isDynamic() {
		return true;
	}

	@Override
	public int getNextPos(ConfIndex confIndex, RCs rcs) {

		int bestPos = -1;
		double bestScore = optimizer.initDouble();

		for (int i=0; i<confIndex.numUndefined; i++) {

			int pos = confIndex.undefinedPos[i];
			double score = scorePos(confIndex, rcs, pos);

			if (optimizer.isBetter(score, bestScore)) {
				bestScore = score;
				bestPos = pos;
			}
		}

		if (bestPos >= 0) {
			return bestPos;
		}

		// sometimes, all the positions have infinite energies
		// so just pick one arbitrarily
		return confIndex.undefinedPos[0];
	}

	double scorePos(ConfIndex confIndex, RCs rcs, int pos) {

		// check all the RCs at this pos and aggregate the energies
		double parentScore = confIndex.node.getScore();
		double reciprocalSum = 0;
		double maxLower = Double.NEGATIVE_INFINITY;
		double minLower = Double.POSITIVE_INFINITY;
		for (int rc : rcs.get(pos)) {
			double childScore = gscorer.calcDifferential(confIndex, rcs, pos, rc)
				+ hscorer.calcDifferential(confIndex, rcs, pos, rc);
			maxLower = Math.max(maxLower, childScore);
			minLower = Math.min(minLower, childScore);
		}

		return minLower -maxLower;
	}
}
