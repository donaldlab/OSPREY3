package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInter;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.ematrix.compiled.EmatCalculator;
import edu.duke.cs.osprey.energy.compiled.*;
import edu.duke.cs.osprey.gpu.Structs;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static edu.duke.cs.osprey.tools.Log.log;

public class TestMARKStarCCS {
    //public static ConfSpace space = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.tiny.complex.ccsx"));
    public static ConfSpace space = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6ov7.small.complex.ccsx"));

    public void runMARKStarBound(ConfSpace confSpace, double epsilon){

        Parallelism parallelism = Parallelism.makeCpu(8);
        TaskExecutor tasks = parallelism.makeTaskExecutor();

        // Define the minimizing energy calculator
        //ConfEnergyCalculator confEcalc = ConfEnergyCalculator.makeBest(confSpace, parallelism);
        ConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);
        ConfEnergyCalculatorAdapter minimizingEcalc = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
                .setPosInterDist(PosInterDist.DesmetEtAl1992)
                .setMinimize(true)
                .setIncludeStaticStatic(false)
                .build();
        ConfEnergyCalculatorAdapter rigidEcalc = new ConfEnergyCalculatorAdapter.Builder(confEcalc, tasks)
                .setPosInterDist(PosInterDist.DesmetEtAl1992)
                .setMinimize(false)
                .setIncludeStaticStatic(false)
                .build();

        EnergyMatrix minEmat = new EmatCalculator.Builder(confEcalc)
                .setPosInterDist(PosInterDist.DesmetEtAl1992)
                .setMinimize(true)
                .setIncludeStaticStatic(false)
                .build()
                .calc(tasks);
        EnergyMatrix rigidEmat = new EmatCalculator.Builder(confEcalc)
                .setPosInterDist(PosInterDist.DesmetEtAl1992)
                .setMinimize(false)
                .setIncludeStaticStatic(false)
                .build()
                .calc(tasks);

        RCs rcs = space.seqSpace().makeWildTypeSequence().makeRCs(space);

        UpdatingEnergyMatrix correctionEmat = new UpdatingEnergyMatrix(space, minEmat);

        MARKStarBoundFastQueues bound = new MARKStarBoundFastQueues(
                space,
                rigidEmat,
                minEmat,
                minimizingEcalc,
                rcs,
                parallelism
        );

        bound.init(epsilon);
        bound.setCorrections(correctionEmat);

        bound.compute();
    }

    @Test
    public void testBound(){
        runMARKStarBound(space, 0.68);
    }

    @Test
    public void testEnergies(){

        log("reading conf space ...");
        var file = new File("/home/graham/IdeaProjects/OSPREY3/src/test/resources/confSpaces/2RL0.A.ccsx");
        var confSpace = ConfSpace.fromBytes(FileTools.readFileBytes(file));
        log("done!");

        // get the wild type conformation, if possible
        var conf = Arrays.stream(confSpace.positions)
                .mapToInt(pos ->
                        Arrays.stream(pos.confs)
                                .filter(c -> c.id.startsWith("wt-"))
                                .mapToInt(c -> c.index)
                                .findFirst()
                                .orElse(0) // otherwise, pick an arbitrary conformation
                )
                .toArray();
        log(Arrays.toString(conf));

        var ecalc = new CPUConfEnergyCalculator(confSpace);
        log("making energy calculator ...");
        log("Ecalc: %s", ecalc.getClass().getSimpleName());

        //var inters = PosInterDist.all(confSpace, new int[]{28,40,28, -1});
        var inters = PosInterDist.all(confSpace, conf);

        var minConf = ecalc.minimize(conf, inters);
        log("\tmin:         %f", minConf.energy);
        log("\tcalc:         %f", ecalc.calcEnergy(minConf.coords, inters));
        log("\n");

        var sumsCalc = new ArrayList<Double>();
        var sumsMin = new ArrayList<Double>();
        for (var inter : inters){
            var interList = new ArrayList<PosInter>();
            interList.add(inter);
            var calc = ecalc.calcEnergy(minConf.coords,interList);
            var min = ecalc.minimizeEnergy(conf,interList);
            log("\tcalcEnergy %f", calc);
            //log("\tminEnergy %f", min);
            sumsMin.add(min);
            sumsCalc.add(calc);

        }
        log("\t sum (Calc) %f", sumsCalc.stream().reduce(0.0, Double::sum));
        log("\t sum (Min) %f", sumsMin.stream().reduce(0.0, Double::sum));


        log("done");
    }

    //workflow:
    // get energiedcoords object
    // use energiedcoords object with Ecalc (compiled) to compute interactions
    @Test
    public void testTriples(){
        log("reading conf space ...");
        var file = new File("/home/graham/IdeaProjects/OSPREY3/src/test/resources/confSpaces/2RL0.A.ccsx");
        var confSpace = ConfSpace.fromBytes(FileTools.readFileBytes(file));
        log("done!");

        // get the wild type conformation, if possible
        var conf = Arrays.stream(confSpace.positions)
                .mapToInt(pos ->
                        Arrays.stream(pos.confs)
                                .filter(c -> c.id.startsWith("wt-"))
                                .mapToInt(c -> c.index)
                                .findFirst()
                                .orElse(0) // otherwise, pick an arbitrary conformation
                )
                .toArray();
        log(Arrays.toString(conf));

        var ecalc = new CPUConfEnergyCalculator(confSpace);
        log("making energy calculator ...");
        log("Ecalc: %s", ecalc.getClass().getSimpleName());

        //var inters = PosInterDist.all(confSpace, new int[]{28,40,28, -1});
        var inters = PosInterDist.all(confSpace, conf);

        var minConf = ecalc.minimize(conf, inters);
        log("\tmin:         %f", minConf.energy);
        log("\tcalc:         %f", ecalc.calcEnergy(conf, inters));
        log("\n");

        var triple = PosInterDist.DesmetEtAl1992.unweightedTripleCorrection(confSpace, null, 2, 28, 1, 40, 0, 28);
        var calc = ecalc.calcEnergy(conf, triple);
        var min = ecalc.calcEnergy(minConf.coords, triple);
        log("\ttriple corr (orig)%f", calc);
        log("\ttriple corr (min)%f", min);
        log("\ttriple corr (diff)%f", min-calc);

        var dubsTest = new ArrayList<PosInter>();
        var dubs = new ArrayList<List<PosInter>>();
        dubs.add(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 2, 28, 1, 40));
        dubsTest.addAll(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 2, 28, 1, 40));
        dubs.add(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 1, 40, 0, 28));
        dubsTest.addAll(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 1, 40, 0, 28));
        dubs.add(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 2, 28, 0, 28));
        dubsTest.addAll(PosInterDist.DesmetEtAl1992.pair(confSpace, null, 2, 28, 0, 28));

        var sumsCorrection = new ArrayList<Double>();
        for (var inter : dubs){
            var orig = ecalc.calcEnergy(conf, inter);
            var min2 = ecalc.calcEnergy(minConf.coords, inter);
            sumsCorrection.add(min2-orig);
            log("\torig %f", orig);
            log("\tmin %f", min2);
            log("\tcorr %f", min2-orig);

        }
        log("\t sum (Correction) %f", sumsCorrection.stream().reduce(0.0, Double::sum));

        log("\t sum (Correction(test)) %f", ecalc.calcEnergy(minConf.coords, dubsTest) - ecalc.calcEnergy(conf, dubsTest));



    }
}
