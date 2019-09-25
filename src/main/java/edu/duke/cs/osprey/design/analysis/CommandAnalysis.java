package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

public interface CommandAnalysis extends PartitionFunction.ConfListener {
    void printResults();
}
