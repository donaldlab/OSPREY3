package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

public interface OptimizableAStarNode {

	double getGScore();
	void setGScore(double val);
	double getHScore();
	void setHScore(double val);
	int getLevel();

	default double getScore() {
		return getGScore() + getHScore();
	}

	default double getGScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getGScore(), optimizer);
	}
	default void setGScore(double val, MathTools.Optimizer optimizer) {
		setGScore(Tools.optimizeScore(val, optimizer));
	}

	default double getHScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getHScore(), optimizer);
	}
	default void setHScore(double val, MathTools.Optimizer optimizer) {
		setHScore(Tools.optimizeScore(val, optimizer));
	}

	default double getScore(MathTools.Optimizer optimizer) {
		return Tools.optimizeScore(getScore(), optimizer);
	}

	public static class Tools {

		public static double optimizeScore(double score, MathTools.Optimizer optimizer) {
			switch (optimizer) {
				case Minimize: return score; // the pq is naturally a min-heap
				case Maximize: return -score; // negate the score so the pq acts like a max-heap
				default: throw new UnpossibleError();
			}
		}
	}
}
