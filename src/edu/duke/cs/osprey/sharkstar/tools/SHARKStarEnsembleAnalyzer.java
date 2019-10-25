package edu.duke.cs.osprey.sharkstar.tools;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.math.BigDecimal;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.sharkstar.EnergyMatrixCorrector.makeTuple;

public class SHARKStarEnsembleAnalyzer {
    Map<RCTuple, BigDecimal> tupleErrors = new ConcurrentHashMap<>();
    private final ConfEnergyCalculator minimizingEcalc;
    private final EnergyMatrix minimizingEmat;
    BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    public SHARKStarEnsembleAnalyzer(ConfEnergyCalculator minimizingEcalc, EnergyMatrix minimizingEmat) {
        this.minimizingEcalc = minimizingEcalc;
        this.minimizingEmat = minimizingEmat;

    }

    public void analyzeFullConf(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf){
        /* What data do we want? */
        computeSubsetEnergies(minimizingEcalc, new RCTuple(conf.getAssignments()),
                analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All), analysis.epmol.energy);

    }

    private void computeSubsetEnergies(ConfEnergyCalculator ecalc, RCTuple tuple, EnergyMatrix fullConfBreakdown,
                                       double fullConfMinEnergy) {
        System.out.println("Analyzing "+tuple.stringListing());
        if(tuple.size() < 2)
            return;
        double energyFromFullConfMinimization = fullConfBreakdown.getInternalEnergy(tuple);
        double fullMinimizedTupleEnergy = ecalc.calcEnergy(tuple).energy;
        double pairwiseMinimizedenergy = minimizingEmat.getInternalEnergy(tuple);
        if(!tupleErrors.containsKey(tuple))
            tupleErrors.put(tuple, BigDecimal.ZERO);
        tupleErrors.put(tuple, tupleErrors.get(tuple)
                .add(bcalc.calc(pairwiseMinimizedenergy).subtract(bcalc.calc(fullMinimizedTupleEnergy))));
        for(int removedPositionIndex = 0; removedPositionIndex < tuple.pos.size(); removedPositionIndex++){
            RCTuple newTuple = tuple.subtractMember(removedPositionIndex);
            computeSubsetEnergies(ecalc, newTuple, fullConfBreakdown, fullConfMinEnergy);
        }
    }

    public void printStats() {
        Map<RCTuple, BigDecimal> sortedTuples = tupleErrors
                .entrySet()
                .stream()
                .sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
                .collect(
                        Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2,
                                LinkedHashMap::new));
        for(RCTuple tuple: sortedTuples.keySet()) {
            System.out.printf("%s:%12.6e\n", tuple, tupleErrors.get(tuple));
        }
    }
}
