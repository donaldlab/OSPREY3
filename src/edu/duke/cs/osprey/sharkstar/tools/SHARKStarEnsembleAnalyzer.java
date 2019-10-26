package edu.duke.cs.osprey.sharkstar.tools;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import org.jetbrains.annotations.NotNull;

import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

public class SHARKStarEnsembleAnalyzer {
    Map<RCTuple, TupleDiff> tupleErrors = new ConcurrentHashMap<>();
    private final ConfEnergyCalculator minimizingEcalc;
    private final EnergyMatrix minimizingEmat;
    BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private int numConfs = 0;
    private BigDecimal ensemblePfuncMinimized = BigDecimal.ZERO;
    private BigDecimal ensemblePfuncLowerBound = BigDecimal.ZERO;

    public SHARKStarEnsembleAnalyzer(ConfEnergyCalculator minimizingEcalc, EnergyMatrix minimizingEmat) {
        this.minimizingEcalc = minimizingEcalc;
        this.minimizingEmat = minimizingEmat;

    }

    public void analyzeFullConf(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf){
        /* What data do we want? */
        numConfs++;
        ensemblePfuncMinimized = ensemblePfuncMinimized.add(bcalc.calc(analysis.epmol.energy));
        ensemblePfuncLowerBound = ensemblePfuncLowerBound.add(bcalc.calc(analysis.score));
        computeSubsetEnergies(minimizingEcalc, new RCTuple(conf.getAssignments()),
                analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All),
                    analysis.epmol.energy, conf.getScore(), new HashSet<RCTuple>());

    }

    private void computeSubsetEnergies(ConfEnergyCalculator ecalc, RCTuple tuple, EnergyMatrix fullConfBreakdown,
                                       double fullMinimizedEnergy, double fullConfPairwiseBound,
                                       Set<RCTuple> processedPartialConfs) {
        if(tuple.size() < 2 || processedPartialConfs.contains(tuple))
            return;
        processedPartialConfs.add(tuple);
        double energyFromFullConfMinimization = fullConfBreakdown.getInternalEnergy(extractPositions(tuple));
        double fullMinimizedTupleEnergy = ecalc.calcEnergy(tuple).energy;
        double pairwiseMinimizedenergy = minimizingEmat.getInternalEnergy(tuple);
        double correctedEnergy = fullMinimizedEnergy + pairwiseMinimizedenergy - energyFromFullConfMinimization;
        if(!tupleErrors.containsKey(tuple))
            tupleErrors.put(tuple, new TupleDiff(tuple));
        tupleErrors.get(tuple).addDiff(fullMinimizedTupleEnergy, pairwiseMinimizedenergy, fullMinimizedEnergy, correctedEnergy);
        for(int removedPositionIndex = 0; removedPositionIndex < tuple.pos.size(); removedPositionIndex++){
            RCTuple newTuple = tuple.subtractMember(removedPositionIndex);
            computeSubsetEnergies(ecalc, newTuple, fullConfBreakdown, fullMinimizedEnergy, fullConfPairwiseBound,
                processedPartialConfs);
        }
    }

    private RCTuple extractPositions(RCTuple tuple) {
        RCTuple newTuple = new RCTuple();
        for(int pos: tuple.pos)
            newTuple = newTuple.addRC(pos, 0);
        return newTuple;
    }

    public void printStats() {
        System.out.printf("Stats for ensemble of %d confs: full pfunc = %12.6e%n", numConfs, ensemblePfuncMinimized);
        Map<RCTuple, TupleDiff> sortedTuples = tupleErrors
                .entrySet()
                .stream()
                .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
                .collect(
                        Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
                                LinkedHashMap::new));
        for(RCTuple tuple: sortedTuples.keySet()) {
            System.out.println(tupleErrors.get(tuple));
        }
    }

    private class TupleDiff implements Comparable<TupleDiff> {
        double totalEnergyDiff = 0;
        int numUniqueConfs = 0;
        BigDecimal totalPfuncDiff = BigDecimal.ZERO;
        public final RCTuple tuple;

        public TupleDiff(RCTuple tuple) {
            this.tuple = tuple;
        }
        @Override
        public int compareTo(@NotNull TupleDiff o) {
            return totalPfuncDiff.compareTo(o.totalPfuncDiff);
        }

        public void addDiff(double fullMinimizedTupleEnergy, double pairwiseMinimizedenergy,
                            double fullMinimizedEnergy, double correctedEnergy) {
            totalEnergyDiff += correctedEnergy - fullMinimizedEnergy;
            numUniqueConfs++;
            totalPfuncDiff = totalPfuncDiff.add(bcalc.calc(correctedEnergy).subtract(bcalc.calc(fullMinimizedEnergy)));
        }

        public String toString() {
            return String.format("%s: (%f) across %d confs -> %12.6e",
                    tuple, totalEnergyDiff, numUniqueConfs, totalPfuncDiff);
        }
    }
}
