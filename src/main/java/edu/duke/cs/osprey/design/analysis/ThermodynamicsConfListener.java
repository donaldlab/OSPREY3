package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.Constants;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.ExpFunction;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

public class ThermodynamicsConfListener implements CommandAnalysis {

    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();
    private BigDecimal pFuncLowerBound;
    private BigDecimal pFuncUpperBound;
    private ExpFunction eFn = new ExpFunction();
    private static BoltzmannCalculator bCalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        confs.add((ConfSearch.EnergiedConf) conf);
    }

    @Override
    public void finished(PartitionFunction pfunc) {
        var values = pfunc.getValues();
        pFuncUpperBound = values.calcUpperBound();
        pFuncLowerBound = pfunc.getValues().calcLowerBound();
    }

    private static BigMath makeBigMath() {
        return new BigMath(PartitionFunction.decimalPrecision);
    }

    private static BigDecimal getBoundedProbability(ConfSearch.EnergiedConf conf, BigDecimal bound) {
        return makeBigMath()
                .set(bCalc.calc(conf.getEnergy()))
                .div(bound)
                .get();
    }

    private BigDecimal getUpperBoundProbability(ConfSearch.EnergiedConf conf) {
        return getBoundedProbability(conf, pFuncLowerBound);
    }

    private BigDecimal getLowerBoundProbability(ConfSearch.EnergiedConf conf) {
        return getBoundedProbability(conf, pFuncUpperBound);
    }

    public BigDecimal getUpperBoundEnthalpy() {
        var math = makeBigMath();
        return confs.stream()
                .map(conf -> math.set(getLowerBoundProbability(conf)).mult(conf.getEnergy()).get())
                .reduce(BigDecimal::add)
                .orElseThrow();
    }

    public BigDecimal getLowerBoundEnthalpy() {
        var upperBoundProbabilitySum = confs.stream()
                .map(this::getUpperBoundProbability)
                .reduce(BigDecimal::add)
                .orElseThrow();

        var restOfProbability = makeBigMath().set(1).sub(upperBoundProbabilitySum);
        var restOfProbabilityTimesWorstEnergyConf = restOfProbability.mult(confs.get(confs.size() - 1).getEnergy()).get();

        return upperBoundProbabilitySum.add(restOfProbabilityTimesWorstEnergyConf);
    }

    public BigDecimal getUpperBoundEntropy() {
        return makeBigMath().set(Constants.R)
                .mult(eFn.log(pFuncUpperBound))
                .add(getUpperBoundEnthalpy().doubleValue() / Constants.T)
                .get();
    }

    public BigDecimal getLowerBoundEntropy() {
        return makeBigMath().set(Constants.R)
                .mult(eFn.log(pFuncLowerBound))
                .add(getLowerBoundEnthalpy().doubleValue() / Constants.T)
                .get();
    }

    @Override
    public void printResults() {
        System.out.println(String.format("Enthalpy[%.04f - %.04f]", getLowerBoundEnthalpy(), getUpperBoundEnthalpy()));
        System.out.println(String.format("Entropy[%.04f - %.04f]", getLowerBoundEntropy(), getUpperBoundEntropy()));
    }
}
