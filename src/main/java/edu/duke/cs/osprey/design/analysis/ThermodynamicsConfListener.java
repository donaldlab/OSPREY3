package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.Constants;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import static ch.obermuhlner.math.big.DefaultBigDecimalMath.log;

public class ThermodynamicsConfListener implements CommandAnalysis {

    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();
    private BigDecimal pFuncLowerBound;
    private BigDecimal pFuncUpperBound;
    private static BoltzmannCalculator bCalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        confs.add((ConfSearch.EnergiedConf) conf);
    }

    /*
    @Override
    public void finished(PartitionFunction pfunc) {
        var values = pfunc.getValues();
        pFuncUpperBound = values.calcUpperBound();
        pFuncLowerBound = pfunc.getValues().calcLowerBound();
    }
     */

    private static BigMath makeBigMath() {
        return new BigMath(PartitionFunction.decimalPrecision);
    }

    private static BigDecimal getBoundedProbability(ConfSearch.EnergiedConf conf, BigDecimal bound) {
        return bCalc.calc(conf.getEnergy()).divide(bound, RoundingMode.HALF_EVEN);
    }

    private ProbabilityEnergyTuple getConfProperties(ConfSearch.EnergiedConf conf) {
        return new ProbabilityEnergyTuple(conf, getLowerBoundProbability(conf), getUpperBoundProbability(conf), BigDecimal.valueOf(conf.getEnergy()));
    }

    private Stream<ProbabilityEnergyTuple> confPropStream() {
        return confs.stream().parallel().map(this::getConfProperties);
    }

    private BigDecimal getUpperBoundProbability(ConfSearch.EnergiedConf conf) {
        return getBoundedProbability(conf, pFuncLowerBound);
    }

    private BigDecimal getLowerBoundProbability(ConfSearch.EnergiedConf conf) {
        return getBoundedProbability(conf, pFuncUpperBound);
    }

    private BigDecimal getUpperBoundEnthalpy() {
        return confPropStream()
                .map(props -> props.lowerBoundProbability().multiply(props.energy()))
                .reduce(BigDecimal::add)
                .orElseThrow();
    }

    private BigDecimal getLowerBoundEnthalpy() {
        var part1 = confPropStream()
                .map(prop -> prop.upperBoundProbability().multiply(prop.energy()))
                .reduce(BigDecimal::add)
                .orElseThrow();

        var upperBoundProbabilitySum = confPropStream()
                .map(ProbabilityEnergyTuple::upperBoundProbability)
                .reduce(BigDecimal::add)
                .orElseThrow();

        var part2 = BigDecimal.ONE
                .subtract(upperBoundProbabilitySum)
                .multiply(confPropStream().reduce((first, second) -> second).orElseThrow().energy());

        return part1.add(part2);
    }

    private BigDecimal getUpperBoundEntropy() {
        return makeBigMath().set(Constants.R)
                .mult(log(pFuncUpperBound))
                .add(getUpperBoundEnthalpy().doubleValue() / Constants.T)
                .get();
    }

    private BigDecimal getLowerBoundEntropy() {
        return makeBigMath().set(Constants.R)
                .mult(log(pFuncLowerBound))
                .add(getLowerBoundEnthalpy().doubleValue() / Constants.T)
                .get();
    }

    @Override
    public void printResults() {
        System.out.println(String.format("Enthalpy[%.04f - %.04f]", getLowerBoundEnthalpy(), getUpperBoundEnthalpy()));
        System.out.println(String.format("Entropy[%.04f - %.04f]", getLowerBoundEntropy(), getUpperBoundEntropy()));
    }
}

record ProbabilityEnergyTuple(ConfSearch.EnergiedConf conf, BigDecimal lowerBoundProbability, BigDecimal upperBoundProbability, BigDecimal energy) {
}
