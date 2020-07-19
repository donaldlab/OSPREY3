package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Thin wrapper class to play nice with BBK* and MSK*
 */
public class SingleSequenceSHARKStarBound implements PartitionFunction {
    public int maxMinimizations = 1;
    private MultiSequenceSHARKStarBound multiSequenceSHARKStarBound;
    public final Sequence sequence;
    private MultiSequenceSHARKStarBound multisequenceBound;
    private Status status;
    private MultiSequenceSHARKStarBound.Values values;
    private int numConfsEvaluated = 0;
    public final BigInteger numConformations;
    public SHARKStarQueue fringeNodes;
    public SHARKStarQueue internalQueue;
    public SHARKStarQueue leafQueue;
    private double sequenceEpsilon = 1;
    private BigDecimal finishedNodeZ = BigDecimal.ZERO;
    public final RCs seqRCs;
    public final State state;

    //debug variable
    public Set<MultiSequenceSHARKStarNode> finishedNodes = new HashSet<>();
    private boolean errors;

    public SingleSequenceSHARKStarBound(MultiSequenceSHARKStarBound multiSequenceSHARKStarBound, Sequence seq, MultiSequenceSHARKStarBound sharkStarBound) {
        this.multiSequenceSHARKStarBound = multiSequenceSHARKStarBound;
        this.sequence = seq;
        this.multisequenceBound = sharkStarBound;
        this.seqRCs = seq.makeRCs(sharkStarBound.confSpace);
        this.numConformations = seqRCs.getNumConformations();
        this.fringeNodes = new SHARKStarQueue(seq);
        this.internalQueue = new SHARKStarQueue(seq);
        this.leafQueue = new SHARKStarQueue(seq);

        this.state = new State();
        this.state.lowerBound = BigDecimal.ZERO;
        this.state.upperBound = BigDecimal.ZERO;
    }

    public void addFinishedNode(MultiSequenceSHARKStarNode node) {
        finishedNodeZ = finishedNodeZ.add(node.getUpperBound(sequence));
        //System.out.println("Adding "+node.toSeqString(sequence)+" to finished set");
        if(finishedNodes.contains(node))
            System.err.println("Dupe node addition.");
        finishedNodes.add(node);
    }

    @Override
    public void setReportProgress(boolean val) {
        multisequenceBound.setReportProgress(val);
    }

    @Override
    public void setConfListener(ConfListener val) {
        multisequenceBound.setConfListener(val);
    }

    @Override
    public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
        init(confSearch, null, numConfsBeforePruning, targetEpsilon);
    }

    @Override
    public void init(ConfSearch upperBoundConfs, ConfSearch lowerBoundConfs, BigInteger numConfsBeforePruning, double targetEpsilon) {
        values = new MultiSequenceSHARKStarBound.Values();
        setStatus(Status.Estimating);
    }


    @Override
    public void setStabilityThreshold(BigDecimal stabilityThreshold) {
        multisequenceBound.setStabilityThreshold(stabilityThreshold);
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
            //multisequenceBound.computeForSequence(maxNumConfs, this);
        multisequenceBound.computeForSequenceParallel(maxNumConfs, this);
        updateBound();
        if (getSequenceEpsilon() < multiSequenceSHARKStarBound.targetEpsilon) {
            setStatus(Status.Estimated);
            if (values.qstar.compareTo(BigDecimal.ZERO) == 0) {
                setStatus(Status.Unstable);
            }
        }
    }

    @Override
    public void compute() {
        compute(Integer.MAX_VALUE);
    }

    @Override
    public Result makeResult() {
        //multiSequenceSHARKStarBound.lowerReduction_ConfUpperBound = multiSequenceSHARKStarBound.rootNode.getLowerBound(sequence)
                //.subtract(multiSequenceSHARKStarBound.startLowerBound).subtract(multiSequenceSHARKStarBound.lowerReduction_FullMin);
        // Calculate the lower bound z reductions from conf upper bounds, since we don't explicitly record these
        //multiSequenceSHARKStarBound.upperReduction_ConfLowerBound = multiSequenceSHARKStarBound.startUpperBound.subtract(multiSequenceSHARKStarBound.rootNode.getUpperBound(sequence))
                //.subtract(multiSequenceSHARKStarBound.upperReduction_FullMin).subtract(multiSequenceSHARKStarBound.upperReduction_PartialMin);

        Result result = new Result(getStatus(), getValues(), getNumConfsEvaluated());
        /*
        result.setWorkInfo(numPartialMinimizations, numConfsScored,minList);
        result.setZInfo(lowerReduction_FullMin, lowerReduction_ConfUpperBound, upperReduction_FullMin, upperReduction_PartialMin, upperReduction_ConfLowerBound);
        result.setOrigBounds(startUpperBound, startLowerBound);
        result.setTimeInfo(stopwatch.getTimeNs());
        result.setMiscInfo(new BigDecimal(rootNode.getNumConfs()));
        */
        return result;
    }

    public BigDecimal getUpperBound(){
        return values.pstar;
    }

    public BigDecimal getLowerBound(){
        return values.qstar;
    }

    public void updateBound() {
        //multiSequenceSHARKStarBound.rootNode.computeEpsilonErrorBounds(sequence);
        BigDecimal lastUpper = values.pstar;
        BigDecimal lastLower = values.qstar;
        BigDecimal upperBound = fringeNodes.getPartitionFunctionUpperBound()
                .add(internalQueue.getPartitionFunctionUpperBound())
                .add(leafQueue.getPartitionFunctionUpperBound())
                .add(finishedNodeZ);
        BigDecimal lowerBound = fringeNodes.getPartitionFunctionLowerBound()
                .add(internalQueue.getPartitionFunctionLowerBound())
                .add(leafQueue.getPartitionFunctionLowerBound())
                .add(finishedNodeZ);
        if(MathTools.isLessThan(lowerBound, lastLower) && ! MathTools.isRelativelySame(lowerBound, lastLower, PartitionFunction.decimalPrecision, 10)) {
            System.err.println("Bounds getting looser. Lower bound is getting lower...");
            errors = true;
        }
        if(MathTools.isGreaterThan(upperBound, lastUpper) && ! MathTools.isRelativelySame(upperBound, lastUpper, PartitionFunction.decimalPrecision, 10)) {
            System.err.println("Bounds getting looser. Upper bound is getting bigger...");
            errors = true;
        }
        values.pstar = upperBound;
        values.qstar = lowerBound;
        values.qprime = upperBound;
        if (upperBound.subtract(lowerBound).compareTo(BigDecimal.ONE) < 1) {
            sequenceEpsilon = 0;
        } else {
            sequenceEpsilon = upperBound.subtract(lowerBound)
                    .divide(upperBound, RoundingMode.HALF_UP).doubleValue();
        }
    }
    List<SeqSpace.ResType> getRTs(SimpleConfSpace.Position confPos, SeqAStarNode.Assignments assignments) {

        // TODO: pre-compute this somehow?
        SeqSpace seqSpace = sequence.seqSpace;

        // map the conf pos to a sequence pos
        SeqSpace.Position seqPos = seqSpace.getPosition(confPos.resNum);
        if (seqPos != null) {

            Integer assignedRT = assignments.getAssignment(seqPos.index);
            if (assignedRT != null) {
                // use just the assigned res type
                return Collections.singletonList(seqPos.resTypes.get(assignedRT));
            } else {
                // use all the res types at the pos
                return seqPos.resTypes;
            }

        } else {

            // immutable position, use all the res types (should just be one)
            assert (confPos.resTypes.size() == 1);

            // use the null value to signal there's no res type here
            return Collections.singletonList(null);
        }
    }

    List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, SeqSpace.ResType rt, SHARKStar.State state) {
        // TODO: pre-compute this somehow?
        if (rt != null) {
            // mutable pos, grab the RCs that match the RT
            return pos.resConfs.stream()
                    .filter(rc -> rc.template.name.equals(rt.name))
                    .collect(Collectors.toList());
        } else {
            // immutable pos, use all the RCs
            return pos.resConfs;
        }
    }

    public boolean nonZeroLower() {
        return this.fringeNodes.getPartitionFunctionLowerBound().compareTo(BigDecimal.ZERO) > 0
                || internalQueue.getPartitionFunctionLowerBound().compareTo(BigDecimal.ZERO) > 0
                || leafQueue.getPartitionFunctionLowerBound().compareTo(BigDecimal.ZERO) > 0 ;
    }

    @Override
    public void printStats() {
        multisequenceBound.printEnsembleAnalysis();
    }
    public boolean errors() {
        return errors;
    }

    public double getSequenceEpsilon() {
        return sequenceEpsilon;
    }


    public void setStatus(Status status) {
        this.status = status;
    }

    public boolean isEmpty() {
        return fringeNodes.isEmpty()
                && internalQueue.isEmpty()
                && leafQueue.isEmpty();
    }

    public static class State{
        BigDecimal upperBound;
        BigDecimal lowerBound;
        long numEnergiedConfs = 0;

        BigDecimal getUpperBound(){
            return upperBound;
        }

        BigDecimal getLowerBound(){
            return lowerBound;
        }

        double calcDelta() {
            BigDecimal upperBound = getUpperBound();
            if (MathTools.isZero(upperBound) || MathTools.isInf(upperBound)) {
                return 1.0;
            }
            return new BigMath(PartitionFunction.decimalPrecision)
                    .set(upperBound)
                    .sub(getLowerBound())
                    .div(upperBound)
                    .get()
                    .doubleValue();
        }
    }
}
