/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests.variationalkstar;

import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import javax.management.RuntimeErrorException;

/**
 *
 * @author hmn5
 */
public class VariationalPartFuncTests {

    static String[] dirNums = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    static double epsilon = 0.1;
    static boolean verbose = true;

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("VariationalKStar/PartitionFunction")) {
            throw new Error("This test was designed to be run in test/VariationalKStar/PartitionFunction folder\n\tcwd: " + path);
        }

//        String[] dirNums = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
//            "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
//            "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

//         String[] dirNums = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10"};
        String command = args[2];
        if (command.equalsIgnoreCase("large")) {
            String subDir = "LargeTest/4HEM/";
            for (String dir : dirNums) {
                String run = subDir + dir + "_Run/";
                ConfigFileParser cfp = setupRun(run);
                SearchProblem searchProb = cfp.getSearchProblem();
                double pruningInt = 50;
                precomputeMatrices(searchProb, cfp, pruningInt);

                PartFuncTree tree = new PartFuncTree(searchProb.emat, searchProb.pruneMat);
                Stopwatch.start();
                double logZ = tree.computeEpsilonApprox(0.1);
                Stopwatch.stop();
                double time = Stopwatch.getTimeMs();
                String fileName = run + "results.txt";
                double epsilonReached = tree.effectiveEpsilon;
                double confSpace = getLogConfSpace(searchProb.pruneMat);
                writeOutResults(fileName, logZ, epsilonReached, time, confSpace);
            }
        } else if (command.equalsIgnoreCase("small")) {
            String run = "SingleTest/4HEM/PartFuncTest/";
            ConfigFileParser cfp = setupRun(run);
            SearchProblem searchProb = cfp.getSearchProblem();
            double pruningInt = 50;
            precomputeMatrices(searchProb, cfp, pruningInt);

            PartFuncTree tree = new PartFuncTree(searchProb.emat, searchProb.pruneMat);
            Stopwatch.start();
            double logZ = tree.computeEpsilonApprox(0.1);
            Stopwatch.stop();
            double time = Stopwatch.getTimeMs();
            String fileName = run + "results.txt";
            double epsilonReached = tree.effectiveEpsilon;
            double confSpace = getLogConfSpace(searchProb.pruneMat);
            printResults(fileName, logZ, epsilonReached, time, confSpace);

            double correctLogZ = 701.3538838283972;
            ExpFunction ef = new ExpFunction();

            double lowerBound = ef.log(ef.exp(correctLogZ).multiply(new BigDecimal(1-epsilon))).doubleValue();
            double upperBound = ef.log(ef.exp(correctLogZ).multiply(new BigDecimal(1+epsilon))).doubleValue();
            if (lowerBound > logZ || upperBound < logZ){
                System.out.println("ERROR: LogZ is not within "+epsilon+" accuracty");
                System.out.println("Lower Bound: "+lowerBound);
                System.out.println("Uppder Bound: "+upperBound);
                System.out.println("LogZ: "+logZ);
            }
            else{
                System.out.println("Lower Bound: "+lowerBound);
                System.out.println("Uppder Bound: "+upperBound);
                System.out.println();
                System.out.println("Test Passed");
            }
        } else {
            throw new RuntimeException("Command should either be large or small, not " + command);
        }
    }

    private static void writeOutResults(String filename, double logZ, double epsilon,
            double time, double confSpace) throws Exception {
        File results = new File(filename);
        try {
            FileWriter fw = new FileWriter(results);
            fw.write("LogConfSpace: " + confSpace + "\n");
            fw.write("Epsilon: " + epsilon + "\n");
            fw.write("LogZ: " + logZ + "\n");
            fw.write("Time: " + time + "\n");
            fw.close();
            if (verbose) {
                System.out.println(filename + " statistics");
                System.out.println("LogConfSpace: " + confSpace);
                System.out.println("Epsilon: " + epsilon);
                System.out.println("LogZ: " + logZ);
                System.out.println("Time: " + time);
                System.out.println();
                System.out.println();
            }
        } catch (IOException e) {
            throw new RuntimeException();
        }
    }

    private static void printResults(String filename, double logZ, double epsilon,
            double time, double confSpace) {
        System.out.println(filename + " statistics");
        System.out.println("LogConfSpace: " + confSpace);
        System.out.println("Epsilon: " + epsilon);
        System.out.println("LogZ: " + logZ);
        System.out.println("Time: " + time);
    }

    private static ConfigFileParser setupRun(String run) {
        String dee = run + "DEE.cfg";
        String sys = run + "System.cfg";
        String kstar = run + "KStar.cfg";

        String[] newArgs = {"-c", kstar, " ", dee, sys};
        ConfigFileParser cfp = new ConfigFileParser(newArgs);
        String pdbFile = cfp.params.getValue("PDBNAME");
        String pathToPDB = run + pdbFile;
        cfp.params.setValue("PDBNAME", pathToPDB);
        String name = cfp.params.getValue("RUNNAME");
        String pathName = run + name;
        cfp.params.setValue("RUNNAME", pathName);
        cfp.loadData();

        return cfp;
    }

    static void precomputeMatrices(SearchProblem searchSpace, ConfigFileParser cfp, double pruningInterval) {
            //Precalculate TupleMatrices needed for GMEC computation.  Some of these may already be computed.  
        //All of these matrices except the basic pairwise energy matrix are pruning-dependent:
        //we can prune conformations whose energies are within pruningInterval
        //of the lowest pairwise lower bound

        //First calculate the pairwise energy matrix, if not already present
        searchSpace.loadEnergyMatrix();

        //Doing competitor pruning now
        //will limit us to a smaller, but effective, set of competitors in all future DEE
        if (searchSpace.competitorPruneMat == null) {
            System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
            PruningControl compPruning = cfp.setupPruning(searchSpace, 0, false, false);
            compPruning.setOnlyGoldstein(true);
            compPruning.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
            System.out.println("COMPETITOR PRUNING DONE");
        }

        //Next, do DEE, which will fill in the pruning matrix
        PruningControl pruning = cfp.setupPruning(searchSpace, pruningInterval, false, false);

        pruning.prune();//pass in DEE options, and run the specified types of DEE            
    }

    private static double getLogConfSpace(PruningMatrix pruneMat) {
        double logConfSpace = 0;
        for (int pos = 0; pos < pruneMat.numPos(); pos++) {
            logConfSpace += Math.log(pruneMat.unprunedRCsAtPos(pos).size());
        }
        return logConfSpace;
    }

}
