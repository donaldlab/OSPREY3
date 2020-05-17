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
import java.util.*;
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
    public SHARKStarQueue_refactor internalQueue; // queue containing the internal node fringe
    public SHARKStarQueue_refactor leafQueue; // queue containing the leaf node fringe
    private double sequenceEpsilon = 1;
    private BigDecimal finishedNodeZ = BigDecimal.ZERO;
    public final RCs seqRCs;
    public State state;

    //debug variable
    public Set<SHARKStarNode> finishedNodes = new HashSet<>();
    private boolean errors;

    public static class State{
        private double lowerBound;
        private double upperBound;
        private BoltzmannCalculator bc;

        public State(BoltzmannCalculator bc){
            this.bc = bc;
            lowerBound = Double.NEGATIVE_INFINITY;
            upperBound = Double.POSITIVE_INFINITY;
        }

        public void setInitialBounds(double lowerBound, double upperBound){
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }

        public void addDiffToLB(double diff){
            this.lowerBound = bc.calc_EDiff(lowerBound, diff);
        }
        public void subtractDiffFromUB(double diff){
            this.upperBound = bc.calc_ESum(upperBound, diff);
        }
    }

    public SingleSequenceSHARKStarBound_refactor(MultiSequenceSHARKStarBound_refactor multiSequenceSHARKStarBound, Sequence seq, BoltzmannCalculator bc) {
        this.multisequenceBound = multiSequenceSHARKStarBound;
        this.sequence = seq;
        this.seqRCs = seq.makeRCs(multisequenceBound.getConfSpace());
        this.numConformations = seqRCs.getNumConformations();
        this.internalQueue = new SHARKStarQueue_refactor(seq, bc);
        this.leafQueue = new SHARKStarQueue_refactor(seq, bc);
        this.state = new State(multiSequenceSHARKStarBound.bc);
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
        //System.out.println("Adding "+node.toString()+" to finished set for " + sequence);
        if(finishedNodes.contains(node))
            System.err.println("Dupe node addition.");
        finishedNodes.add(node);
    }

    public double calcEpsilon(){
        return getSequenceEpsilon();
        /*
        BigDecimal ZUB = calcZBound(e -> e.getFreeEnergyLB(sequence));
        BigDecimal ZLB = calcZBound(e -> e.getFreeEnergyUB(sequence));
        return ZUB.subtract(ZLB, decimalPrecision)
                .divide(ZUB, decimalPrecision).doubleValue();

         */
    }

    public double calcEpsilonDirect(){
        BigDecimal ZUB = calcZBound(e -> e.getFreeEnergyLB(sequence));
        BigDecimal ZLB = calcZBound(e -> e.getFreeEnergyUB(sequence));
        return ZUB.subtract(ZLB, decimalPrecision)
                .divide(ZUB, decimalPrecision).doubleValue();
    }

    public BigDecimal calcZUBDirect(){
        return Stream.of(internalQueue, leafQueue, finishedNodes)
                .flatMap(Collection::stream)
                .parallel()
                .map(e -> multisequenceBound.bc.calc(e.getFreeEnergyLB(sequence)))
                .reduce(BigDecimal.ZERO,BigDecimal::add);
    }

    public BigDecimal calcZLBDirect(){
        return Stream.of(internalQueue, leafQueue, finishedNodes)
                .flatMap(Collection::stream)
                .parallel()
                .map(e -> multisequenceBound.bc.calc(e.getFreeEnergyUB(sequence)))
                .reduce(BigDecimal.ZERO,BigDecimal::add);
    }

    public double calcEBound(Function<SHARKStarNode, Double> energyMapper){
        // If the queues are empty, treat this as positive infinity
        if (internalQueue.isEmpty() && leafQueue.isEmpty() && finishedNodes.isEmpty()){
            return Double.POSITIVE_INFINITY;
        }else {
            Optional<Double> minElement = Stream.of(internalQueue, leafQueue, finishedNodes)
                    .flatMap(Collection::stream)
                    .parallel()
                    .map(energyMapper)
                    .min(Double::compareTo);

            if (minElement.isPresent()) {
                Stream<Double> allNodesEnergyLB = Stream.of(internalQueue, leafQueue, finishedNodes)
                        .flatMap(Collection::stream)
                        .parallel()
                        .map(energyMapper);

                if (!allNodesEnergyLB.isParallel()) {
                    System.err.println("ERROR: stream is not parallel");
                }

                return multisequenceBound.bc.logSumExpStream(allNodesEnergyLB, minElement.get());
            }
            else{
                throw new RuntimeException("Could not find min element of stream");
            }
        }

    }

    public BigDecimal calcZBound(Function <SHARKStarNode, Double> energyMapper){
        return multisequenceBound.bc.calc(calcEBound(energyMapper));
    }

    public double getELB(){
        return state.lowerBound;
    }

    public double getEUB(){
        return state.upperBound;
    }

    public BigDecimal getZLB(){
        return multisequenceBound.bc.calc(getEUB());
    }
    public BigDecimal getZUB(){
        return multisequenceBound.bc.calc(getELB());
    }

    @Deprecated
    public double getSequenceEpsilon(){
        BigDecimal ZUB = state.bc.calc(state.lowerBound);
        BigDecimal ZLB = state.bc.calc(state.upperBound);
        if (ZUB.subtract(ZLB).compareTo(BigDecimal.ONE) < 1) {
            return 0;
        } else {
            return ZUB.subtract(ZLB)
                    .divide(ZUB, RoundingMode.HALF_UP).doubleValue();
        }
    }

    public boolean isEmpty() {
        return internalQueue.isEmpty()
                && leafQueue.isEmpty();
    }

}
