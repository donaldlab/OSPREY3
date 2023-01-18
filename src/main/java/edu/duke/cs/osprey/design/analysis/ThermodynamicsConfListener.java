package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.Constants;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.auto.value.AutoValue;

import static ch.obermuhlner.math.big.DefaultBigDecimalMath.log;

public class ThermodynamicsConfListener implements CommandAnalysis {

    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();
    private BigDecimal pFuncLowerBound;
    private BigDecimal pFuncUpperBound;
    private static BoltzmannCalculator bCalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private int maxNumConfs;

    public ThermodynamicsConfListener() { this(-1); }

    public ThermodynamicsConfListener(int maxConfs) {
        maxNumConfs = maxConfs;
    }

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        if (maxNumConfs > 0 && confs.size() > maxNumConfs) {
            return; // don't save any further conformations
        }

        confs.add((ConfSearch.EnergiedConf) conf);
    }

    public void setPFuncBounds(BigDecimal lowerBound, BigDecimal upperBound) {
        this.pFuncLowerBound = lowerBound;
        this.pFuncUpperBound = upperBound;
    }

    private static BigMath makeBigMath() {
        return new BigMath(PartitionFunction.decimalPrecision);
    }

    private static BigDecimal getBoundedProbability(ConfSearch.EnergiedConf conf, BigDecimal bound) {
        return bCalc.calc(conf.getEnergy()).divide(bound, RoundingMode.HALF_EVEN);
    }

    private ProbabilityEnergyTuple getConfProperties(ConfSearch.EnergiedConf conf) {
        return ProbabilityEnergyTuple.create(conf, getLowerBoundProbability(conf), getUpperBoundProbability(conf));
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
        // write out the conformations and energies of the ensemble (if we've got some limit)
        if (maxNumConfs > 0) {

            System.out.println("energy,lprob,uprob,conf");
            confPropStream().forEach(t -> {
                var prettyRcs = Arrays.stream(t.conf().getAssignments())
                        .mapToObj(String::valueOf)
                        .collect(Collectors.joining(";"));

                BigDecimal lProb = t.lowerBoundProbability();
                BigDecimal uProb = t.upperBoundProbability();
                BigDecimal energy = t.energy();

                System.out.printf("%f,%f,%f,%s%n", energy, lProb, uProb, prettyRcs);
            });
        }

        System.out.printf("Enthalpy[%.04f - %.04f]%n", getLowerBoundEnthalpy(), getUpperBoundEnthalpy());
        System.out.printf("Entropy[%.04f - %.04f]%n", getLowerBoundEntropy(), getUpperBoundEntropy());
    }
}

@AutoValue
abstract class ProbabilityEnergyTuple {

    static ProbabilityEnergyTuple create(ConfSearch.EnergiedConf conf, BigDecimal lowerBoundProbability, BigDecimal upperBoundProbability) {
        return new AutoValue_ProbabilityEnergyTuple(conf, lowerBoundProbability, upperBoundProbability, BigDecimal.valueOf(conf.getEnergy()));
    }

    abstract ConfSearch.EnergiedConf conf();
    abstract BigDecimal lowerBoundProbability();
    abstract BigDecimal upperBoundProbability();
    abstract BigDecimal energy();
}
