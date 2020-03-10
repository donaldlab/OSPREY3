package edu.duke.cs.osprey.design.analysis;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class EnsembleAnalysisListener implements CommandAnalysis {

    private final int maxConfs;
    private final PartitionFunction pfunc;
    private final ConfEnergyCalculator eCalc;
    private final String outputDir;
    private List<ConfSearch.ScoredConf> confs = new ArrayList<>();

    public EnsembleAnalysisListener(PartitionFunction pfunc, ConfEnergyCalculator eCalc, int maxConfs, String outputDir) {
        this.pfunc = pfunc;
        this.eCalc = eCalc;
        this.maxConfs = maxConfs;
        this.outputDir = Paths.get(outputDir).toAbsolutePath().toString();
    }

    @Override
    public void printResults() {
        var path = String.format("%s%s%s", outputDir, File.separator, "conf-*.pdb");
        System.out.println(String.format("Writing ensemble out to %s", path));
        var analysis = new ConfAnalyzer(eCalc).analyzeEnsemble(confs.iterator(), maxConfs);
        analysis.writePdbs(path);
    }

    @Override
    public void onConf(ConfSearch.ScoredConf conf) {
        if (confs.size() > maxConfs) {
            return;
        }
        confs.add(conf);
    }

    @Override
    public void finished(PartitionFunction pfunc) {

    }
}
