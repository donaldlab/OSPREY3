package edu.duke.cs.osprey.astar.conf.scoring;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public interface ScorerFactory {
    AStarScorer make(EnergyMatrix emat);
}
