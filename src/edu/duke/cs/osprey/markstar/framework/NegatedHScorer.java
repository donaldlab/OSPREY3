package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;

public class NegatedHScorer implements AStarScorer {
    private AStarScorer sourceScorer;
    public NegatedHScorer(AStarScorer scorer) {
        this.sourceScorer = scorer;
    }

    @Override
    public AStarScorer make() {
        return new NegatedHScorer(sourceScorer);
    }

    @Override
    public double calc(ConfIndex confIndex, RCs rcs) {
        return -sourceScorer.calc(confIndex, rcs);
    }
}
