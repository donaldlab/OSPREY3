package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.Collection;
import java.util.HashSet;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Thin wrapper class to play nice with BBK* and MSK*
 *
 * TODO: Determine whether I need to go back to SHARKStarQueues for speed optimization purposes
 */
public class SingleSequenceSHARKStarBound_refactor implements PartitionFunction {
    public int maxMinimizations = 1;
    private MultiSequenceSHARKStarBound_refactor multiSequenceSHARKStarBound;
    public final Sequence sequence;
    private MultiSequenceSHARKStarBound_refactor multisequenceBound;
    private Status status;
    private MultiSequenceSHARKStarBound_refactor.Values values;
    private int numConfsEvaluated = 0;
    public final BigInteger numConformations;
    public SHARKStarQueue_refactor fringeNodes;
    public SHARKStarQueue_refactor internalQueue;
    public SHARKStarQueue_refactor leafQueue;
    private double sequenceEpsilon = 1;
    private BigDecimal finishedNodeZ = BigDecimal.ZERO;
    public final RCs seqRCs;

    //debug variable
    private Set<SHARKStarNode> finishedNodes = new HashSet<>();
    private boolean errors;

    public SingleSequenceSHARKStarBound_refactor(MultiSequenceSHARKStarBound_refactor multiSequenceSHARKStarBound, Sequence seq, BoltzmannCalculator bc) {
        this.multisequenceBound = multiSequenceSHARKStarBound;
        this.sequence = seq;
        this.seqRCs = seq.makeRCs(multisequenceBound.getConfSpace());
        this.numConformations = seqRCs.getNumConformations();
        this.fringeNodes = new SHARKStarQueue_refactor(seq, bc);
        this.internalQueue = new SHARKStarQueue_refactor(seq, bc);
        this.leafQueue = new SHARKStarQueue_refactor(seq, bc);
    }
    @Override
    public void setReportProgress(boolean val) {
        multisequenceBound.setReportProgress(val);
    }

    @Override
    public void setConfListener(ConfListener val) {
        multisequenceBound.setConfListener(val);
    }

    public void setStatus(Status val){
        this.status = val;
    }

    @Override
    public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
        values = new MultiSequenceSHARKStarBound.Values();
        setStatus(Status.Estimating);
    }

    @Override
    public Status getStatus() {
        return this.status;
    }

    @Override
    public Values getValues() {
        return this.values;
    }

    @Override
    public int getParallelism() {
        return multisequenceBound.getParallelism();
    }

    @Override
    public int getNumConfsEvaluated() {
        return numConfsEvaluated;
    }

    @Override
    public void compute(int maxNumConfs) {
        multisequenceBound.computeForSequence(maxNumConfs, this);
    }

    public void compute(){
        compute(Integer.MAX_VALUE);
    }

    public void addFinishedNode(SHARKStarNode node) {
        System.out.println("Adding "+node.toString()+" to finished set for " + sequence);
        if(finishedNodes.contains(node))
            System.err.println("Dupe node addition.");
        finishedNodes.add(node);
    }

    public double calcEpsilon(){
        BigDecimal ZUB = calcZBound(e -> e.getFreeEnergyLB(sequence));
        return ZUB.subtract(calcZBound(e -> e.getFreeEnergyUB(sequence)), decimalPrecision)
                .divide(ZUB, decimalPrecision).doubleValue();

    }

    public BigDecimal calcZBound(Function<SHARKStarNode, Double> energyMapper){
        Stream<Double> allNodesEnergyLB = Stream.of(fringeNodes, internalQueue, leafQueue, finishedNodes)
                .flatMap(Collection::stream)
                .parallel()
                .map(energyMapper);

       if (!allNodesEnergyLB.isParallel()){
            System.err.println("ERROR: stream is not parallel");
        }

        return multisequenceBound.bc.calc(multisequenceBound.bc.logSumExpStream(allNodesEnergyLB));

    }

    public double getSequenceEpsilon(){
        BigDecimal upperBound = fringeNodes.getPartitionFunctionUpperBound()
                .add(internalQueue.getPartitionFunctionUpperBound())
                .add(leafQueue.getPartitionFunctionUpperBound())
                .add(finishedNodeZ);
        BigDecimal lowerBound = fringeNodes.getPartitionFunctionLowerBound()
                .add(internalQueue.getPartitionFunctionLowerBound())
                .add(leafQueue.getPartitionFunctionLowerBound())
                .add(finishedNodeZ);

        if (upperBound.subtract(lowerBound).compareTo(BigDecimal.ONE) < 1) {
            return 0;
        } else {
            return upperBound.subtract(lowerBound)
                    .divide(upperBound, RoundingMode.HALF_UP).doubleValue();
        }
    }

    public boolean isEmpty() {
        return fringeNodes.isEmpty()
                && internalQueue.isEmpty()
                && leafQueue.isEmpty();
    }

}
