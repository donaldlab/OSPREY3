/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dataAnalysis.variationalpartfunc;

import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.DiscretePartFunc;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author hmn5
 */
public class DataLargeRuns {

    static String[] dirNums4HEM = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    static String[] dirNums4LAJ = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26"};

    static String[] dirNums3GXU = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    static double[] epsilons = {0.8, 0.6, 0.4, 0.2, 0.1};
    static double maxTime = 3600000;
    static boolean verbose = true;
    static boolean printStats = false;

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
        if (!path.endsWith("DataAnalysis/VariationalPartFunc/LargeRuns")) {
            throw new Error("This test was designed to be run in DataAnalysis/VariationalPartFunc/LargeRuns folder\n\tcwd: " + path);
        }
        if (args.length >= 2 && args[0].equalsIgnoreCase("KSTAR")) {
            System.out.println("Performing KStar Analysis");
            if (args[1].equalsIgnoreCase("4HEM")) {
                String[] subDirsToRun = {"4HEM/Lovell/"};
                runAnalysisKStar(subDirsToRun);
                runAnalysisKStarUnbound(subDirsToRun);
            } else if (args[1].equalsIgnoreCase("4LAJ")) {
                String[] subDirsToRun = {"4LAJ/Lovell/"};
                runAnalysisKStar(subDirsToRun);
                runAnalysisKStarUnbound(subDirsToRun);
            } else if (args[1].equalsIgnoreCase("3GXU")) {
                String[] subDirsToRun = {"3GXU/Lovell/"};
                runAnalysisKStar(subDirsToRun);
                runAnalysisKStarUnbound(subDirsToRun);
            }
        } else {
            runAnalysisVariationalPartFunc();
            runAnalysisVariationalPartFuncUnbound();
        }
    }

    private static void runAnalysisVariationalPartFunc() throws Exception {
        PartFuncTree.verbose = false;
        String[] subDirs = {"4LAJ/Lovell/", "4HEM/Lovell/", "3GXU/Lovell/"};
        String[] runList;
        for (String subDir : subDirs) {
            if (subDir.contains("4LAJ")) {
                runList = dirNums4LAJ;
            } else if (subDir.contains("4HEM")) {
                runList = dirNums4HEM;
            } else {
                runList = dirNums3GXU;
            }
            for (String dir : runList) {
                String run = subDir + dir + "_Run/";
                System.out.println("Running Analysis in: " + run);
                ConfigFileParser cfp = setupRun(run);
                SearchProblem searchProb = cfp.getSearchProblem();
                double pruningInt = 50;
                precomputeMatrices(searchProb, cfp, pruningInt);
                for (double epsilon : epsilons) {
                    System.out.println("\t Epsilon: " + epsilon);
                    PartFuncTree tree = new PartFuncTree(searchProb.emat, searchProb.pruneMat);
                    Stopwatch.start();
                    double logZ = tree.computeEpsilonApprox(epsilon, maxTime);
                    Stopwatch.stop();
                    double time = Stopwatch.getTimeMs();
                    boolean didFinish = tree.timeOut;
                    double effectiveEpsilon = tree.effectiveEpsilon;
                    String fileName = run + "Data/varPF_" + epsilon + "_" + 1 + ".txt";
                    double epsilonReached = tree.effectiveEpsilon;
                    double confSpace = getLogConfSpace(searchProb.pruneMat);
                    writeOutResults(fileName, logZ, epsilonReached, time, confSpace, didFinish, effectiveEpsilon);
                }
                System.out.println();
            }
        }
    }

    private static void runAnalysisVariationalPartFuncUnbound() throws Exception {
        PartFuncTree.verbose = false;
        String[] subDirs = {"4LAJ/Lovell/", "4HEM/Lovell/", "3GXU/Lovell/"};
        String[] runList;
        for (String subDir : subDirs) {
            if (subDir.contains("4LAJ")) {
                runList = dirNums4LAJ;
            } else if (subDir.contains("4HEM")) {
                runList = dirNums4HEM;
            } else {
                runList = dirNums3GXU;
            }
            for (String dir : runList) {
                String run = subDir + dir + "_Run/";
                System.out.println("Running Analysis in: " + run);
                ConfigFileParser cfp = setupRun(run);
                SearchProblem[] searchProbs = cfp.getMSDSearchProblems();
                for (int i = 1; i < 3; i++) {
                    SearchProblem searchProb = searchProbs[i];
                    double pruningInt = 50;
                    precomputeMatrices(searchProb, cfp, pruningInt);
                    for (double epsilon : epsilons) {
                        System.out.println("\t Epsilon: " + epsilon);
                        PartFuncTree tree = new PartFuncTree(searchProb.emat, searchProb.pruneMat);
                        Stopwatch.start();
                        double logZ = tree.computeEpsilonApprox(epsilon, maxTime);
                        Stopwatch.stop();
                        double time = Stopwatch.getTimeMs();
                        boolean didFinish = tree.timeOut;
                        double effectiveEpsilon = tree.effectiveEpsilon;
                        String fileName = run + "Data_Unbound_" + (i - 1) + "/varPF_" + epsilon + "_" + 1 + ".txt";
                        double epsilonReached = tree.effectiveEpsilon;
                        double confSpace = getLogConfSpace(searchProb.pruneMat);
                        writeOutResults(fileName, logZ, epsilonReached, time, confSpace, didFinish, effectiveEpsilon);
                    }
                    System.out.println();
                }
            }
        }
    }

    private static void runAnalysisKStar(String[] subDirToRun) throws Exception {
        PartFuncTree.verbose = false;
        String[] runList;
        for (String subDir : subDirToRun) {
            if (subDir.contains("4LAJ")) {
                runList = dirNums4LAJ;
            } else if (subDir.contains("4HEM")) {
                runList = dirNums4HEM;
            } else {
                runList = dirNums3GXU;
            }
            for (String dir : runList) {
                String run = subDir + dir + "_Run/";
                System.out.println("Running KStar Analysis in: " + run);
                ConfigFileParser cfp = setupRun(run);
                SearchProblem searchProb = cfp.getSearchProblem();
                double pruningInt = 50;
                precomputeMatrices(searchProb, cfp, pruningInt);
                for (double epsilon : epsilons) {
                    System.out.println("\t Epsilon: " + epsilon);
                    DiscretePartFunc dpf = new DiscretePartFunc(searchProb.emat, searchProb.pruneMat, epsilon, maxTime);
                    String fileName = run + "Data/kstarPF_" + epsilon + "_" + 1 + ".txt";
                    double confSpace = getLogConfSpace(searchProb.pruneMat);
                    double logZ = dpf.getLogZ();
                    double time = dpf.totalTime;
                    boolean didFinish = dpf.finishedInTime;
                    double effectivEpsilon = dpf.effectiveEpsilonReached;
                    writeOutResults(fileName, logZ, epsilon, time, confSpace, didFinish, effectivEpsilon);
                }
                System.out.println();
            }
        }
    }

    private static void runAnalysisKStarUnbound(String[] subDirToRun) throws Exception {
        PartFuncTree.verbose = false;
        String[] runList;
        for (String subDir : subDirToRun) {
            if (subDir.contains("4LAJ")) {
                runList = dirNums4LAJ;
            } else if (subDir.contains("4HEM")) {
                runList = dirNums4HEM;
            } else {
                runList = dirNums3GXU;
            }
            for (String dir : runList) {
                String run = subDir + dir + "_Run/";
                System.out.println("Running KStar Analysis in: " + run);
                ConfigFileParser cfp = setupRun(run);
                SearchProblem[] searchProbs = cfp.getMSDSearchProblems();
                for (int i = 1; i < 3; i++) {
                    SearchProblem searchProb = searchProbs[i];
                    double pruningInt = 50;
                    precomputeMatrices(searchProb, cfp, pruningInt);
                    for (double epsilon : epsilons) {
                        System.out.println("\t Epsilon: " + epsilon);
                        DiscretePartFunc dpf = new DiscretePartFunc(searchProb.emat, searchProb.pruneMat, epsilon, maxTime);
                        String fileName = run + "Data_Unbound_" + (i - 1) + "/kstarPF_" + epsilon + "_" + 1 + ".txt";
                        double confSpace = getLogConfSpace(searchProb.pruneMat);
                        double logZ = dpf.getLogZ();
                        double time = dpf.totalTime;
                        boolean didFinish = dpf.finishedInTime;
                        double effectivEpsilon = dpf.effectiveEpsilonReached;
                        writeOutResults(fileName, logZ, epsilon, time, confSpace, didFinish, effectivEpsilon);
                    }
                    System.out.println();
                }
            }
        }
    }

    private static void writeOutResults(String filename, double logZ, double epsilon,
            double time, double confSpace, boolean didFinish, double effectiveEpsilon) throws Exception {
        File results = new File(filename);
        try {
            FileWriter fw = new FileWriter(results);
            fw.write("LogConfSpace: " + confSpace + "\n");
            fw.write("Epsilon: " + epsilon + "\n");
            fw.write("LogZ: " + logZ + "\n");
            fw.write("Time: " + time + "\n");
            if (didFinish) {
                fw.write("Finished: true");
            } else {
                fw.write("Finished: false");
            }
            fw.write("EffectiveEpsilon: " + effectiveEpsilon);
            fw.close();
            if (printStats) {
                printResults(filename, logZ, epsilon, time, confSpace);
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
        System.out.println();
        System.out.println();
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
            PruningControl compPruning = cfp.setupPruning(searchSpace, 0, false, false);
            compPruning.setOnlyGoldstein(true);
            compPruning.prune();
            searchSpace.competitorPruneMat = searchSpace.pruneMat;
            searchSpace.pruneMat = null;
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
