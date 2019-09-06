package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.SeqSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.math.BigDecimal;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

public class SHARKSeqHScorer implements SeqAStarScorer {
    // TODO: make configurable
    private static final int upperBatchSize = 1000;
    private static final int numLowerBatches = 1;
    public final RTs rts = null;
    private final SeqSpace seqSpace;
    private final Iterable<SHARKStar.WeightedState> objectives;
    private final SimpleConfSpace fullConfSpace;
    private final EnergyMatrix rigidEnergyMatrix;
    private final EnergyMatrix minimizedEnergyMatrix;
    BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private BigDecimal flexibleUpper;
    private BigDecimal flexibleLower;

    public SHARKSeqHScorer(SeqSpace sequenceSpace, SimpleConfSpace fullConfSpace,
                           EnergyMatrix rigidEnergyMatrix, EnergyMatrix minimizedEnergyMatrix) {
        this(sequenceSpace, fullConfSpace, rigidEnergyMatrix, minimizedEnergyMatrix, null);
    }

    public SHARKSeqHScorer(SeqSpace sequenceSpace, SimpleConfSpace fullConfSpace,
                           EnergyMatrix rigidEnergyMatrix,
                           EnergyMatrix minimizedEnergyMatrix,
                           Iterable<SHARKStar.WeightedState> objectives) {
        this.objectives = objectives;
        this.seqSpace = sequenceSpace;
        this.fullConfSpace = fullConfSpace;
        this.rigidEnergyMatrix = rigidEnergyMatrix;
        this.minimizedEnergyMatrix = minimizedEnergyMatrix;
    }

    public SHARKSeqHScorer(SeqSpace sequenceSpace, Iterable<SHARKStar.WeightedState> objectives) {
        this(sequenceSpace, null, null, null, objectives);
    }

    @Override
    public double calc(SeqAStarNode.Assignments assignments) {
        if(objectives != null)
            return calc(objectives, assignments);
        return Double.NaN;
    }

    public double calc(SimpleConfSpace fullConfSpace, SeqAStarNode.Assignments assignments) {
        double lowerBound = 0.0;
        MathTools.Optimizer confOpt = null;
        MathTools.Optimizer seqOpt = null;
        List<SimpleConfSpace.Position> positions = fullConfSpace.positions;
        confOpt = MathTools.Optimizer.Minimize;
        seqOpt = MathTools.Optimizer.Maximize;

        BigDecimal stateSum = BigDecimal.ONE;
        for(SimpleConfSpace.Position pos1 : positions) {
            BigDecimal residueSum = seqOpt.initBigDecimal();
            for (SeqSpace.ResType rt1 : getRTs(pos1, assignments)) {
                BigDecimal AASum = BigDecimal.ZERO;
                for (SimpleConfSpace.ResidueConf rc1 : getRCs(pos1, rt1)) {
                    double rotamerSum = minimizedEnergyMatrix.getEnergy(pos1.index, rc1.index);
                    for(SimpleConfSpace.Position pos2 : positions) {
                        if (pos2.index >= pos1.index)
                            continue;
                        double bestPair = confOpt.initDouble();
                        for (SeqSpace.ResType rt2 : getRTs(pos2, assignments)) {
                            for (SimpleConfSpace.ResidueConf rc2 : getRCs(pos2, rt2)) {
                                bestPair = confOpt.opt(bestPair, minimizedEnergyMatrix.getEnergy(pos1.index, rc1.index,
                                        pos2.index, rc2.index));
                            }
                        }
                        rotamerSum += bestPair;
                    }
                    AASum = AASum.add(bcalc.calc(rotamerSum));
                }
                residueSum = seqOpt.opt(residueSum, AASum);
            }
            stateSum = stateSum.multiply(residueSum);
        }
        System.out.println("Sum for "+assignments+":"+bcalc.freeEnergy(stateSum));
        // add the bound to our estimate of the objective LMFE
        lowerBound += bcalc.freeEnergy(stateSum);

        System.out.println("Lowerbound updated to "+lowerBound);
        return lowerBound;

    }

    public double calc(Iterable<SHARKStar.WeightedState> objective, SeqAStarNode.Assignments assignments) {
        double lowerBound = 0.0;
        MathTools.Optimizer confOpt = null;
        MathTools.Optimizer seqOpt = null;
        for (SHARKStar.WeightedState wstate : objective) {
            List<SimpleConfSpace.Position> positions = wstate.state.confSpace.positions;
            if(wstate.weight < 0) {
                confOpt = MathTools.Optimizer.Maximize;
                seqOpt = MathTools.Optimizer.Minimize;
            }
            if(wstate.weight > 0) {
                confOpt = MathTools.Optimizer.Minimize;
                seqOpt = MathTools.Optimizer.Maximize;
            }

            BigDecimal stateSum = BigDecimal.ONE;
            for(SimpleConfSpace.Position pos1 : wstate.state.confSpace.positions) {
                BigDecimal residueSum = seqOpt.initBigDecimal();
                for (SeqSpace.ResType rt1 : getRTs(pos1, assignments)) {
                    BigDecimal AASum = BigDecimal.ZERO;
                    for (SimpleConfSpace.ResidueConf rc1 : getRCs(pos1, rt1)) {
                        double rotamerSum = wstate.getSingleEnergy(pos1.index, rc1.index);
                        for(SimpleConfSpace.Position pos2 : wstate.state.confSpace.positions) {
                            if (pos2.index >= pos1.index)
                                continue;
                            double bestPair = confOpt.initDouble();
                            for (SeqSpace.ResType rt2 : getRTs(pos2, assignments)) {
                                for (SimpleConfSpace.ResidueConf rc2 : getRCs(pos2, rt2)) {
                                    bestPair = confOpt.opt(bestPair, wstate.getPairEnergy(pos1.index, rc1.index,
                                            pos2.index, rc2.index));
                                }
                            }
                            rotamerSum += bestPair;
                        }
                        AASum = AASum.add(bcalc.calc(rotamerSum));
                    }
                    residueSum = seqOpt.opt(residueSum, AASum);
                }
                stateSum = stateSum.multiply(residueSum);
            }
            System.out.println("Sum for "+assignments+" in state "+wstate.state+":"+bcalc.freeEnergy(stateSum));
            if (wstate.weight > 0) {

                // add the bound to our estimate of the objective LMFE
                lowerBound += wstate.weight * bcalc.freeEnergy(stateSum);

            } else {

                // compute a lower bound on the pfunc value of the multi-sequence conf space
                // add the weighted bound to our estimate of the objective LMFE
                lowerBound += wstate.weight * bcalc.freeEnergy(stateSum);
            }
            System.out.println("Lowerbound updated to "+lowerBound);
        }
        return lowerBound;

    }

    List<SeqSpace.ResType> getRTs(SimpleConfSpace.Position confPos, SeqAStarNode.Assignments assignments) {

        // TODO: pre-compute this somehow?

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

    List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, SeqSpace.ResType rt) {
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
}
