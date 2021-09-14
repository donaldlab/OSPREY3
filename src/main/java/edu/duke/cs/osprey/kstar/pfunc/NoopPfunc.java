package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.parallelism.TaskExecutor;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.time.Duration;

public class NoopPfunc implements PartitionFunction {
    @Override
    public void setReportProgress(boolean val) {
    }

    @Override
    public void setConfListener(ConfListener val) {
    }

    @Override
    public void init(double targetEpsilon) {
    }

    @Override
    public Status getStatus() {
                            return Status.Estimated;
        }

    @Override
    public Values getValues() {
        Values values = new Values();
        values.qstar = BigDecimal.ONE;
        return values;
    }

    @Override
    public int getParallelism() {
                              return 0;
                                       }

    @Override
    public int getNumConfsEvaluated() {
                                    return 0;
                                             }

    @Override
    public void compute(int maxNumConfs) {

    }

    @Override
    public void setStabilityThreshold(BigDecimal stabilityThreshold) {
    }

}
