package edu.duke.cs.osprey.astar;

import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

/**
 * An interface to allow nodes to not implement the standard getScore method, but still
 * make use of the compatible parts of the framework around ConfAStarNode
 */
public interface PartialOptimizableAStarNode {
    double getGScore();

    void setGScore(double val);

    int getLevel();

    default double getGScore(MathTools.Optimizer optimizer) {
        return OptimizableAStarNode.Tools.optimizeScore(getGScore(), optimizer);
    }

    default double getRigidGScore() {
        return getGScore();
    }

    class Tools {

        public static double optimizeScore(double score, MathTools.Optimizer optimizer) {
            switch (optimizer) {
                case Minimize: return score; // the pq is naturally a min-heap
                case Maximize: return -score; // negate the score so the pq acts like a max-heap
                default: throw new UnpossibleError();
            }
        }
    }
}
