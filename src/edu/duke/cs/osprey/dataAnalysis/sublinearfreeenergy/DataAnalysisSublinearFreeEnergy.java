/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.dataAnalysis.sublinearfreeenergy;

import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.SequenceFreeEnergy;
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
public class DataAnalysisSublinearFreeEnergy {
    static String[] dirNums4HEM = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15"};// "16", "17", "18", "19", "20",
//        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    static String[] dirNums4LAJ = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10"};//,
//        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
//        "21", "22", "23", "24", "25", "26"};

    static String[] dirNums3GXU = {"01", "02", "03", "04", "05", "06", "07", "08", "09", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"};
//        "21", "22", "23", "24", "25", "26", "27", "28", "29", "30"};

    static boolean normalizeFreeEnergy = false;

    public static void main(String[] args)
            throws Exception {

        String path = new File("").getAbsolutePath();
//        if (!path.endsWith("DataAnalysis/SublinearFreeEnergy/LargeRuns")) {
//            throw new Error("This test was designed to be run in DataAnalysis/SublinearFreeEnergy/LargeRuns folder\n\tcwd: " + path);
//        }
        boolean doExhaustive = false;
        if (args.length >= 1 && args[0].equalsIgnoreCase("exhaustive")) {
            doExhaustive = true;
        }
        runAnalysisSublinearFreeEnergy(doExhaustive);
    }

    private static void runAnalysisSublinearFreeEnergy(boolean doExhaustive) throws Exception {
        PartFuncTree.verbose = false;
//        String[] subDirs = {"4LAJ/Lovell/", "4HEM/Lovell/", "3GXU/Lovell/"};
//        String[] subDirs = {"3GXU/Lovell/"};
        String[] subDirs = {"4HEM/Lovell/"};
//        String[] subDirs = {"4LAJ/Lovell/"};
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
                if (normalizeFreeEnergy) {
                    SequenceFreeEnergy.normalizeFreeEnergy = true;
                }
                Stopwatch.start();
                SequenceFreeEnergy feTree = new SequenceFreeEnergy(cfp);
                if (doExhaustive) {
                    feTree.exhaustiveSearch();
                } else {
                    feTree.nextConf();
                }
                Stopwatch.stop();
                String[] bestSequence = feTree.bestSequence;
                double time = Stopwatch.getTimeMs();
                double confSpace = getLogConfSpace(searchProb.pruneMat);
                double sequenceSpace = feTree.getLogSequenceSpaceSize();
                if (doExhaustive) {
                    String filename;
                    if (!normalizeFreeEnergy) {
                        filename = run + "Data/exhaustive.txt";
                    } else {
                        filename = run + "Data/exhaustive_normalized.txt";
                    }
                    writeOutResultsExhaustive(filename, bestSequence, time, confSpace, sequenceSpace);
                } else {
                    String filename;
                    if (!normalizeFreeEnergy) {
                        filename = run + "Data/dataSublinear.txt";
                    }
                    else{
                        filename = run + "Data/dataSublinear_normalized.txt";
                    }
                    int numLeafNodesVisited = feTree.numLeafNodesVisited;
                    int numLeafNodesExpanded = feTree.numLeafNodesExpanded;
                    writeOutResults(filename, bestSequence, time, confSpace, sequenceSpace, numLeafNodesVisited, numLeafNodesExpanded);
                }
                System.out.println();
            }
        }
    }

    private static void writeOutResults(String filename, String[] bestSequence,
            double time, double confSpace, double sequenceSpace, int numLeafVis, int numLeafExp) throws Exception {
        File results = new File(filename);
        try {
            FileWriter fw = new FileWriter(results);
            fw.write("LogConfSpace: " + confSpace + "\n");
            fw.write("LogSeqSpace: " + sequenceSpace + "\n");
            fw.write("Time: " + time + "\n");
            fw.write("LeafVisited: " + numLeafVis + "\n");
            fw.write("LeafExpanded: " + numLeafExp + "\n");
            fw.write("Sequence:");
            for (String aa : bestSequence) {
                fw.write(" " + aa);
            }
            fw.close();
        } catch (IOException e) {
            throw new RuntimeException();
        }
    }

    private static void writeOutResultsExhaustive(String filename, String[] bestSequence,
            double time, double confSpace, double sequenceSpace) throws Exception {
        File results = new File(filename);
        try {
            FileWriter fw = new FileWriter(results);
            fw.write("LogConfSpace: " + confSpace + "\n");
            fw.write("LogSeqSpace: " + sequenceSpace + "\n");
            fw.write("Time: " + time + "\n");
            fw.write("Sequence:");
            for (String aa : bestSequence) {
                fw.write(" " + aa);
            }
            fw.close();
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
