package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.TestKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.util.List;

import static edu.duke.cs.osprey.sharkstar.TestSHARKStar.loadFromCFS;
import static edu.duke.cs.osprey.sharkstar.TestSHARKStarBound.makeMultiSequenceSHARKStarPfuncForConfSpace;

public class TestParallelMinimization {

    public static ConfAnalyzer confAnalyzer;
    public static ConfEnergyCalculator minimizingConfEcalc;
    public static EnergyMatrixCorrector energyMatrixCorrector;
    public static ConfSearch.ScoredConf conf;
    public static TestKStar.ConfSpaces confSpaces;

    @BeforeClass
    public static void before() {
        try {
            confSpaces = loadFromCFS("test-resources/3ma2_A_6res_3.157E+06.cfs");
            conf = new ConfSearch.ScoredConf(new int[]{29,4,7,3,4,3}, -61.174);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    @Test
    public void testSingleVsQuad() {
        System.out.println("Testest!");
        Stopwatch times = new Stopwatch().start();
        minimizeOneConformation(1);
        reportTime(times, "Time for one minimization: tt");
        minimizeOneConformation(4);
        reportTime(times, "Time for minimization with four cores: tt");
        computeCorrectionsForOneConformation(1);
        reportTime(times, "Time for corrections with one core: tt");
        computeCorrectionsForOneConformation(4);
        reportTime(times, "Time for corrections with four cores: tt");
    }

    private void reportTime(Stopwatch times, String reportString) {
        times.stop();
        System.out.println(reportString.replace("tt", times.getTime()));
        times.reset();
        times.start();
    }

    private void computeCorrectionsForOneConformation(int numCores) {
        SimpleConfSpace confSpace = confSpaces.complex;
        MultiSequenceSHARKStarBound mssharkbound =
                (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(confSpace,
                        new RCs(confSpace), 0.999, null, numCores);
        energyMatrixCorrector = new EnergyMatrixCorrector(mssharkbound);
        minimizingConfEcalc = mssharkbound.minimizingEcalc;
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
        energyMatrixCorrector.computeEnergyCorrection(analysis, conf, 0);
    }

    private void minimizeOneConformation(int numCores) {
        boolean original = false;
        if(original) {
            minimizingConfEcalc.calcEnergy(conf);
        }
        else {
            SimpleConfSpace confSpace = confSpaces.complex;
            MultiSequenceSHARKStarBound mssharkbound =
                    (MultiSequenceSHARKStarBound) makeMultiSequenceSHARKStarPfuncForConfSpace(confSpace,
                            new RCs(confSpace), 0.999, null, numCores);
            minimizingConfEcalc = mssharkbound.minimizingEcalc;
            confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
            ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
        }
    }
}
