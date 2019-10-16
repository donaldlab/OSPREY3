package edu.duke.cs.osprey.sharkstar.tools;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.math.BigDecimal;
import java.util.HashMap;
import java.util.Map;

import static edu.duke.cs.osprey.sharkstar.EnergyMatrixCorrector.makeTuple;

public class SHARKStarEnsembleAnalyzer {
    Map<RCTuple, BigDecimal> tupleErrors = new HashMap<>();
    ConfEnergyCalculator minimizingEcalc = null;
    EnergyMatrix minimizingEmat = null;
    BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    public SHARKStarEnsembleAnalyzer(ConfEnergyCalculator minimizingEcalc, EnergyMatrix minimizingEmat) {

    }

    public void analyzeFullConf(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                 ConfEnergyCalculator ecalc, double epsilonBound) {
        /* What data do we want? */
        computeSubsetEnergies(ecalc, makeTuple(conf));

    }

    private void computeSubsetEnergies(ConfEnergyCalculator ecalc, RCTuple tuple) {
        double fullMinimizedTupleEnergy = ecalc.calcEnergy(tuple).energy;
        double pairwiseMinimizedenergy = minimizingEmat.getInternalEnergy(tuple);
        if(!tupleErrors.containsKey(tuple))
            tupleErrors.put(tuple, BigDecimal.ZERO);
        tupleErrors.put(tuple, tupleErrors.get(tuple)
                .add(bcalc.calc(fullMinimizedTupleEnergy-pairwiseMinimizedenergy)));
        for(int removedPosition = tuple.size(); removedPosition > 0; removedPosition--){
            computeSubsetEnergies(ecalc, tuple.subtractMember(removedPosition));
        }
    }
}
