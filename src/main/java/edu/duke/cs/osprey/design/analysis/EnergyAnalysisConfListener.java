package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class EnergyAnalysisConfListener implements CommandAnalysis {

    private final ConfEnergyCalculator eCalc;
    private final List<Long> indices;
    private List<ConfSearch.EnergiedConf> confs = new ArrayList<>();
    private long currentIdx;

    public EnergyAnalysisConfListener(ConfEnergyCalculator eCalc, List<Long> indicesToAnalyze) {
        this.eCalc = eCalc;
        this.indices = indicesToAnalyze;
    }

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        if (indices.contains(currentIdx)) {
            confs.add((ConfSearch.EnergiedConf) conf);
            indices.remove(currentIdx);
        }

        currentIdx++;
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
        System.out.println("Energy analysis follows:\n");
        for (var idx : IntStream.range(0, indices.size()).toArray()) {
            var analysis = analyzeConf(idx);
            System.out.println(String.format("Energy Analysis of %d sequence", idx + 1));
            System.out.println(analysis + "\n");

            System.out.println(String.format("Forcefield breakdown: \n%s\n", analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All)));
            System.out.println(String.format("Electrostatic breakdown: \n%s\n", analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.Electrostatics)));
            System.out.println(String.format("van der Waals breakdown: \n%s\n", analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.VanDerWaals)));
            System.out.println(String.format("Solvation breakdown: \n%s\n", analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.Solvation)));
            System.out.println(String.format("Offsets breakdown: \n%s\n", analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.Offsets)));
        }
    }
}
