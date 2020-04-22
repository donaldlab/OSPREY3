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

import static edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound.debug;

class SHARKStarQueue_refactor extends PriorityQueue<SHARKStarNode> {
    private BigDecimal partitionFunctionUpperSum = BigDecimal.ZERO;
    private BigDecimal partitionFunctionLowerSum = BigDecimal.ZERO;
    private final Sequence seq;
    private final BoltzmannCalculator bc;

    public SHARKStarQueue_refactor(Sequence seq, BoltzmannCalculator bc) {
        //TODO: Is this extra logging / exping worth it? Is it slow?
        super((o1, o2) -> -Double.compare(
                bc.calc_lnZDiff(o1.getFreeEnergyLB(seq), o1.getFreeEnergyUB(seq)),
                bc.calc_lnZDiff(o2.getFreeEnergyLB(seq), o2.getFreeEnergyUB(seq))
        ));
        this.seq = seq;
        this.bc = bc;
    }

    public BigDecimal getPartitionFunctionUpperBound() {
        return partitionFunctionUpperSum;
    }

    public BigDecimal getPartitionFunctionLowerBound() {
        return partitionFunctionLowerSum;
    }

    /*
    @Override
    public boolean add(SHARKStarNode node) {
        debugCheck();
        partitionFunctionUpperSum = partitionFunctionUpperSum.add(node.getUpperBound(seq));
        partitionFunctionLowerSum = partitionFunctionLowerSum.add(node.getLowerBound(seq));
        debugCheck();
        return super.add(node);
    }

    @Override
    public SHARKStarNode poll() {
        SHARKStarNode node = super.poll();
        debugCheck();
        partitionFunctionUpperSum = partitionFunctionUpperSum.subtract(node.getUpperBound(seq));
        partitionFunctionLowerSum = partitionFunctionLowerSum.subtract(node.getLowerBound(seq));
        debugCheck();
        return node;
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
        if (!isEmpty() && peek().getLowerBound(seq).compareTo(partitionFunctionLowerSum) > 0)
            System.err.println("The top element is bigger than the entire lower bound sum.");
        assert (sumDifference.compareTo(BigDecimal.ZERO) > 0 || sumDifference.compareTo(BigDecimal.valueOf(1e-5)) <= 0);
        assert (partitionFunctionLowerSum.compareTo(BigDecimal.ZERO) >= 0);
        System.out.println("Queue: bounds " + toString());
        List<MultiSequenceSHARKStarNode> nodes = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            if (isEmpty())
                break;
            MultiSequenceSHARKStarNode next = super.poll();
            System.out.println(next.toSeqString(seq));
            nodes.add(next);
        }
        for (MultiSequenceSHARKStarNode node : nodes)
            super.add(node);
    }

    private void debugCheck() {
        debugCheck(false);
    }
     */
}
