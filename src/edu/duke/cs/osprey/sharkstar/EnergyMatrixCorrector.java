package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupE;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.BatchCorrectionMinimizer;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;

public class EnergyMatrixCorrector {
    private final MultiSequenceSHARKStarBound multiSequenceSHARKStarBound;
    private ConfEnergyCalculator confEcalc;
    private BatchCorrectionMinimizer batcher;

    public EnergyMatrixCorrector(MultiSequenceSHARKStarBound multiSequenceSHARKStarBound) {
        this.multiSequenceSHARKStarBound = multiSequenceSHARKStarBound;
        this.confEcalc = multiSequenceSHARKStarBound.minimizingEcalc;
        this.batcher = new BatchCorrectionMinimizer(confEcalc, multiSequenceSHARKStarBound.correctionMatrix,
                multiSequenceSHARKStarBound.minimizingEmat);
    }

    void computeEnergyCorrection(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                 double epsilonBound) {
        if (conf.getAssignments().length < 3)
            return;
        //System.out.println("Analysis:"+analysis);
        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(multiSequenceSHARKStarBound.getMinimizingEmat());
        Stopwatch correctionTime = new Stopwatch().start();
        //System.out.println("Energy Analysis: "+energyAnalysis);
        //System.out.println("Score Analysis: "+scoreAnalysis);
        EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);
        //System.out.println("Difference Analysis " + diff);
        List<TupE> sortedPairwiseTerms2 = new ArrayList<TupE>();
        for (int pos = 0; pos < diff.getNumPos(); pos++) {
            for (int rc = 0; rc < diff.getNumConfAtPos(pos); rc++) {
                for (int pos2 = 0; pos2 < diff.getNumPos(); pos2++) {
                    for (int rc2 = 0; rc2 < diff.getNumConfAtPos(pos2); rc2++) {
                        if (pos >= pos2)
                            continue;
                        double sum = 0;
                        sum += diff.getOneBody(pos, rc);
                        sum += diff.getPairwise(pos, rc, pos2, rc2);
                        sum += diff.getOneBody(pos2, rc2);
                        TupE tupe = new TupE(new RCTuple(pos, rc, pos2, rc2), sum);
                        sortedPairwiseTerms2.add(tupe);
                    }
                }
            }
        }
        Collections.sort(sortedPairwiseTerms2);

        double threshhold = 0.1;
        double minDifference = 0.9;
        double triplethreshhold = 0.3;
        // storePartialConfCorrections(conf, epsilonBound, diff, sortedPairwiseTerms2, threshhold, minDifference, triplethreshhold);
        storePartialConfCorrections(conf, epsilonBound, diff, sortedPairwiseTerms2, threshhold, minDifference, triplethreshhold,
                4);
        batcher.submit();
        confEcalc.tasks.waitForFinish();
        correctionTime.stop();
    }

    private void storePartialConfCorrections(ConfSearch.ScoredConf conf, double epsilonBound, EnergyMatrix diff,
                                             List<TupE> sortedPairwiseTerms2, double threshhold, double minDifference,
                                             double minTupleDiff, int maxTupleSize) {
        double maxDiff = sortedPairwiseTerms2.get(0).E;
        for (int i = 0; i < sortedPairwiseTerms2.size(); i++) {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if (pairDiff < minDifference && maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            recursePartialCorrection(conf, epsilonBound, diff, minTupleDiff, maxTupleSize, makeTuple(conf, pos1, pos2));
        }
    }

    private void recursePartialCorrection(ConfSearch.ScoredConf conf, double epsilonBound, EnergyMatrix diff, double minTupleDiff, int maxTupleSize, RCTuple curTuple) {
        if(curTuple.size() > maxTupleSize)
            return;
        for (int nextPos = 0; nextPos < diff.getNumPos(); nextPos++) {
            if (curTuple.pos.contains(nextPos))
                continue;
            RCTuple tuple = curTuple.addRC(nextPos, conf.getAssignments()[nextPos]);
            if(tuple.pos.size() > 2 && !multiSequenceSHARKStarBound.correctionMatrix.hasHigherOrderTermFor(tuple)) {
                double tupleBounds = multiSequenceSHARKStarBound.getRigidEmat().getInternalEnergy(tuple) - multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(tuple);
                if (tupleBounds < minTupleDiff)
                    continue;
                multiSequenceSHARKStarBound.minList.set(tuple.size() - 1, multiSequenceSHARKStarBound.minList.get(tuple.size() - 1) + 1);
                batcher.getBatch().addTuple(tuple);
                batcher.submitIfFull();
                multiSequenceSHARKStarBound.setNumPartialMinimizations(multiSequenceSHARKStarBound.getNumPartialMinimizations() + 1);
                multiSequenceSHARKStarBound.getProgress().reportPartialMinimization(1, epsilonBound);
            }
            recursePartialCorrection(conf, epsilonBound, diff, minTupleDiff, maxTupleSize, tuple);
        }
    }

    private void storePartialConfCorrections(ConfSearch.ScoredConf conf, double epsilonBound, EnergyMatrix diff,
                                             List<TupE> sortedPairwiseTerms2, double threshhold, double minDifference,
                                             double triplethreshhold) {
        double maxDiff = sortedPairwiseTerms2.get(0).E;
        for (int i = 0; i < sortedPairwiseTerms2.size(); i++) {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if (pairDiff < minDifference && maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            int localMinimizations = 0;
            for (int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double tupleBounds = multiSequenceSHARKStarBound.getRigidEmat().getInternalEnergy(tuple) - multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(tuple);
                if (tupleBounds < triplethreshhold)
                    continue;
                multiSequenceSHARKStarBound.minList.set(tuple.size() - 1, multiSequenceSHARKStarBound.minList.get(tuple.size() - 1) + 1);
                computeDifference(tuple, multiSequenceSHARKStarBound.getMinimizingEcalc());
                localMinimizations++;
            }
            multiSequenceSHARKStarBound.setNumPartialMinimizations(multiSequenceSHARKStarBound.getNumPartialMinimizations() + localMinimizations);
            multiSequenceSHARKStarBound.getProgress().reportPartialMinimization(localMinimizations, epsilonBound);
        }
    }

    void computeDifference(RCTuple tuple, ConfEnergyCalculator minimizingEcalc) {
        multiSequenceSHARKStarBound.setComputedCorrections(true);
        if (multiSequenceSHARKStarBound.getCorrectedTuples().contains(tuple.stringListing()))
            return;
        multiSequenceSHARKStarBound.getCorrectedTuples().add(tuple.stringListing());
        if (multiSequenceSHARKStarBound.getCorrectionMatrix().hasHigherOrderTermFor(tuple))
            return;
        double tripleEnergy = minimizingEcalc.calcEnergy(tuple).energy;

        double lowerbound = multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(tuple);
        if (tripleEnergy - lowerbound > 0) {
            double correction = tripleEnergy - lowerbound;
            multiSequenceSHARKStarBound.getCorrectionMatrix().setHigherOrder(tuple, correction);
        } else
            System.err.println("Negative correction for " + tuple.stringListing());

    }

    public static RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for (int pos : positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    void processPreminimization(SingleSequenceSHARKStarBound bound, ConfEnergyCalculator ecalc) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        int maxMinimizations = 1;//parallelism.numThreads;
        List<MultiSequenceSHARKStarNode> topConfs = getTopConfs(queue, maxMinimizations);
        // Need at least two confs to do any partial preminimization
        if (topConfs.size() < 2) {
            queue.addAll(topConfs);
            return;
        }
        RCTuple lowestBoundTuple = topConfs.get(0).toTuple();
        RCTuple overlap = findLargestOverlap(lowestBoundTuple, topConfs, 3);
        //Only continue if we have something to minimize
        for (MultiSequenceSHARKStarNode conf : topConfs) {
            RCTuple confTuple = conf.toTuple();
            if (multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(confTuple) == multiSequenceSHARKStarBound.getRigidEmat().getInternalEnergy(confTuple))
                continue;
            multiSequenceSHARKStarBound.setNumPartialMinimizations(multiSequenceSHARKStarBound.getNumPartialMinimizations() + 1);
            multiSequenceSHARKStarBound.minList.set(confTuple.size() - 1, multiSequenceSHARKStarBound.minList.get(confTuple.size() - 1) + 1);
            if (confTuple.size() > 2 && confTuple.size() < RCs.getNumPos()) {
                multiSequenceSHARKStarBound.getMinimizingEcalc().tasks.submit(() -> {
                    computeTupleCorrection(multiSequenceSHARKStarBound.getMinimizingEcalc(), conf.toTuple(), bound.getSequenceEpsilon());
                    return null;
                }, (econf) -> {
                });
            }
        }
        //minimizingEcalc.tasks.waitForFinish();
        if (overlap.size() > 3 && !multiSequenceSHARKStarBound.getCorrectionMatrix().hasHigherOrderTermFor(overlap)
                && multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(overlap) != multiSequenceSHARKStarBound.getRigidEmat().getInternalEnergy(overlap)) {
            multiSequenceSHARKStarBound.getMinimizingEcalc().tasks.submit(() -> {
                computeTupleCorrection(ecalc, overlap, bound.getSequenceEpsilon());
                return null;
            }, (econf) -> {
            });
        }
        queue.addAll(topConfs);
    }

    void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap, double epsilonBound) {
        if (multiSequenceSHARKStarBound.getCorrectionMatrix().hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = multiSequenceSHARKStarBound.getMinimizingEmat().getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        multiSequenceSHARKStarBound.getProgress().reportPartialMinimization(1, epsilonBound);
        if (partiallyMinimizedLower > pairwiseLower)
            synchronized (multiSequenceSHARKStarBound.getCorrectionMatrix()) {
                multiSequenceSHARKStarBound.getCorrectionMatrix().setHigherOrder(overlap, partiallyMinimizedLower - pairwiseLower);
            }
        multiSequenceSHARKStarBound.getProgress().reportPartialMinimization(1, epsilonBound);
    }

    List<MultiSequenceSHARKStarNode> getTopConfs(PriorityQueue<MultiSequenceSHARKStarNode> queue, int numConfs) {
        List<MultiSequenceSHARKStarNode> topConfs = new ArrayList<MultiSequenceSHARKStarNode>();
        while (topConfs.size() < numConfs && !queue.isEmpty()) {
            MultiSequenceSHARKStarNode nextLowestConf = queue.poll();
            topConfs.add(nextLowestConf);
        }
        return topConfs;
    }

    RCTuple findLargestOverlap(RCTuple conf, List<MultiSequenceSHARKStarNode> otherConfs, int minResidues) {
        RCTuple overlap = conf;
        for (MultiSequenceSHARKStarNode other : otherConfs) {
            overlap = overlap.intersect(other.toTuple());
            if (overlap.size() < minResidues)
                break;
        }
        return overlap;

    }
}