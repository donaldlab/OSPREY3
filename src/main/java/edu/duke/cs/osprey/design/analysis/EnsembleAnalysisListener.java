package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.io.File;
import java.nio.file.Paths;

public class EnsembleAnalysisListener implements CommandAnalysis {

    private final ConfEnergyCalculator eCalc;
    private final String outputDir;
    private final EnergiedConfQueue confs;

    public EnsembleAnalysisListener(ConfEnergyCalculator eCalc, int maxConfs, String outputDir) {
        this.eCalc = eCalc;
        this.outputDir = Paths.get(outputDir).toAbsolutePath().toString();
        this.confs = new EnergiedConfQueue(maxConfs);
    }

    @Override
    public void printResults() {
        var path = String.format("%s%s%s", outputDir, File.separator, "conf-*.pdb");
        System.out.println(String.format("Writing ensemble out to %s", path));
        var lst = confs.toOrderedList();
        var analysis = new ConfAnalyzer(eCalc).analyzeEnsemble(lst.iterator(), lst.size());
        analysis.writePdbs(path);
    }

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        var eConf = ((ConfSearch.EnergiedConf) conf);
        confs.add(eConf);
    }
}
