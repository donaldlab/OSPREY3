package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.control.Main;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.util.ArrayList;
import java.util.List;

public class EnergyAnalysisConfListener implements CommandAnalysis {

    private final ConfEnergyCalculator eCalc;
    private final List<Integer> indices;
    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();

    public EnergyAnalysisConfListener(ConfEnergyCalculator eCalc, List<Integer> indicesToAnalyze) {
        this.eCalc = eCalc;
        this.indices = indicesToAnalyze;
    }

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        confs.add((ConfSearch.EnergiedConf) conf);
    }

    @Override
    public void finished(PartitionFunction pfunc) {

    }

    ConfAnalyzer.ConfAnalysis analyzeConf(int idx) {
        final var analyzer = new ConfAnalyzer(eCalc);
        return analyzer.analyze(confs.get(idx));
    }

    @Override
    public void printResults() {

    }
}
