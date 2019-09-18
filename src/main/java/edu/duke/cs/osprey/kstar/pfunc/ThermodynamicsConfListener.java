package edu.duke.cs.osprey.kstar.pfunc;

import edu.duke.cs.osprey.Constants;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.ExpFunction;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

public class ThermodynamicsConfListener implements PartitionFunction.ConfListener {

    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();
    private BigDecimal pFuncValue;
    private ExpFunction eFn = new ExpFunction();
    private BoltzmannCalculator bCalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        confs.add((ConfSearch.EnergiedConf) conf);
    }

    public void setPartitionFunctionValue(BigDecimal value) {
        this.pFuncValue = value;
    }

    private static BigMath makeBigMath() {
        return new BigMath(PartitionFunction.decimalPrecision);
    }


    private BigDecimal getProbability(ConfSearch.EnergiedConf conf) {
        return makeBigMath()
                .set(bCalc.calc(conf.getEnergy()))
                .div(pFuncValue)
                .get();
    }

    public BigDecimal getEnthalpy() {
        var math = makeBigMath();

        return confs.stream()
                .map(conf -> math.set(getProbability(conf)).mult(conf.getEnergy()).get())
                .reduce(BigDecimal::add)
                .orElseThrow();
    }

    public BigDecimal getEntropy() {
        var math = makeBigMath();

        return math.set(Constants.R)
                .mult(eFn.log(pFuncValue))
                .add(getEnthalpy().doubleValue() / Constants.T)
                .get();
    }
}
