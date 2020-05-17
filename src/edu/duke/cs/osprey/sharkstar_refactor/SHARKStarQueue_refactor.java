package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.concurrent.PriorityBlockingQueue;

import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.debug;

class SHARKStarQueue_refactor extends PriorityBlockingQueue<SHARKStarNode> {
    /**
     * TODO: Try to batch the boltzmann calculator calls?
     */
    private BigDecimal partitionFunctionUpperSum = BigDecimal.ZERO;
    private BigDecimal partitionFunctionLowerSum = BigDecimal.ZERO;
    private final Sequence seq;
    private final BoltzmannCalculator bc;

    public SHARKStarQueue_refactor(Sequence seq, BoltzmannCalculator bc) {
        super(100,(o1, o2) -> -Double.compare(o1.getScore(seq), o2.getScore(seq)));
        this.seq = seq;
        this.bc = bc;
    }

    public BigDecimal getPartitionFunctionUpperBound() {
        return partitionFunctionUpperSum;
    }

    public BigDecimal getPartitionFunctionLowerBound() {
        return partitionFunctionLowerSum;
    }

    public String toString() {
        return "" + size() + " nodes->" + new MathTools.BigDecimalBounds(partitionFunctionLowerSum, partitionFunctionUpperSum);

    }

    public void debugCheck(boolean force) {
        if (!debug || !force)
            return;
        BigDecimal sumDifference = partitionFunctionUpperSum.subtract(partitionFunctionLowerSum);
        if (sumDifference.compareTo(BigDecimal.ZERO) < 0 && sumDifference.compareTo(BigDecimal.valueOf(1e-5)) > 0)
            System.err.println("Invalid bounds. Lower bound is greater than upper bound.");
        if (partitionFunctionLowerSum.compareTo(BigDecimal.ZERO) < 0)
            System.err.println("Invalid bounds. Lower bound is less than zero.");
        if (!isEmpty() && peek().getFreeEnergyUB(seq) < bc.freeEnergy(partitionFunctionLowerSum))
            System.err.println("The top element is bigger than the entire lower bound sum.");
        assert (sumDifference.compareTo(BigDecimal.ZERO) > 0 || sumDifference.compareTo(BigDecimal.valueOf(1e-5)) <= 0);
        assert (partitionFunctionLowerSum.compareTo(BigDecimal.ZERO) >= 0);
        System.out.println("Queue: bounds " + toString());
        List<SHARKStarNode> nodes = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            if (isEmpty())
                break;
            SHARKStarNode next = super.poll();
            //this is a good place for more debugging info
            //System.out.println(next.toSeqString(seq));
            nodes.add(next);
        }
        for (SHARKStarNode node : nodes)
            super.add(node);
    }

    private void debugCheck() {
        debugCheck(false);
    }

}
