/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.comets.UpdatedPruningMatrix;
import edu.duke.cs.osprey.partitionfunctionbounds.GumbelMapTree;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.DiscretePartFunc;
import edu.duke.cs.osprey.partitionfunctionbounds.MinSpanningTree;
import edu.duke.cs.osprey.pruning.Pruner;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * @author hmn5
 */
public class VariationalPartFunc {

    ConfigFileParser cfp;
    SearchProblem sp;
    final static double constRT = PoissonBoltzmannEnergy.constRT;
    ExpFunction ef = new ExpFunction();
    boolean testSCMF = true;
    double epsilon;

    public VariationalPartFunc(ConfigFileParser aCFP) {
        this.cfp = aCFP;
        this.cfp.params.setValue("STERICTHRESH", "1000");
//        SearchProblem[] spList = cfp.getMSDSearchProblems();

  /*      for (SearchProblem searchProb : spList) {
            loadEMatandPrune(searchProb, Double.POSITIVE_INFINITY);
        }*/

//        sp = spList[0];
        sp = cfp.getSearchProblem();
        loadEMatandPrune(sp, Double.POSITIVE_INFINITY);
        if (testSCMF) {
//            testSCMF(sp);

//            DiscretePartFunc dpf = new DiscretePartFunc(sp, 0.1);
//            double logZKStar = dpf.getLogZ();
//            System.out.println("KStar LogZ: " + logZKStar);
            prune(sp, 30);
            System.out.println("Finished Pruning");
            /*if (true) {
             partFuncTree tree = new partFuncTree(sp.emat, upm);
             double logZ = tree.computeEpsilonApprox(0.1);
             System.out.println("Epsilon Approx BB: " + logZ);
             }
             */
            if (false) {
                double average = 0.0;
                int averageNodesExpanded = 0;
                int numSamples = 10;
                for (int i = 0; i < numSamples; i++) {
                    GumbelMapTree tree = new GumbelMapTree(sp);
                    tree.nextConf();
                    double score = tree.currentBestFeasibleScore;
                    average += score;
                    double logZestimate = -average / (this.constRT * (i + 1));
                    System.out.println("Number of Nodes Expanded: " + tree.numExpanded);
                    averageNodesExpanded += tree.numExpanded;
                    System.out.println("Current Sample: " + -score / this.constRT);
                    System.out.println("Current Average: " + logZestimate);
                    System.out.println("Current Average Nodes Exp: " + averageNodesExpanded / (i + 1));
                }
                System.out.println("Average Nodes Expanded: " + averageNodesExpanded / numSamples);
                double logZ = -average / (this.constRT * numSamples);
                System.out.println("Gumbel logZ: " + logZ);
                System.exit(0);
            }
            if (true) {
                epsilon = 0.3;
                PartFuncTree tree = new PartFuncTree(sp.emat, sp.pruneMat);
                long startTime = System.currentTimeMillis();
                double maxTime = 3600000;

                System.out.println("Starting Part Func Calculation with Epsilon: " + epsilon);
                double logZ = tree.computeEpsilonApprox(epsilon, maxTime);
                long totalTime = (System.currentTimeMillis() - startTime);
                System.out.println("New Alg Took: " + totalTime + " milliseconds");
//                DiscretePartFunc dfp = new DiscretePartFunc(sp.emat, sp.pruneMat, 0.1, maxTime);
//                System.out.println("DFP: " + dfp.getLogZ());
                String filename = "data_";
                filename += epsilon + "_1Hour.txt";
                File statistics = new File(filename);
                try {
                    FileWriter fw = new FileWriter(statistics);
                    fw.write("LogConfSpace: " + getLogConfSpace(sp.pruneMat) + "\n");
                    fw.write("NewAlgorithm: " + totalTime + "\n");
                    if (tree.timeOut) {
                        fw.write("NewAlgorithm: finished false" + "\n");
                        fw.write("NewAlgorithm: effectiveEpsilon " + tree.effectiveEpsilon + "\n");
                    } else {
                        fw.write("NewAlgorithm: finished true" + "\n");
                    }
                    fw.write("NewAlgorithm: logZ " + logZ + "\n");
                    DiscretePartFunc dfp = new DiscretePartFunc(sp.emat, sp.pruneMat, 0.1, maxTime);
                    System.out.println("DFP: " + dfp.getLogZ());
                    if (dfp.finishedInTime) {
                        fw.write("KStar: finished true" + "\n");
                        fw.write("KStar: totalTime " + dfp.totalTime);
                    } else {
                        fw.write("KStar: finished false" + "\n");
                        fw.write("KStar: effectiveEpsilon " + dfp.effectiveEpsilonReached + "\n");
                        fw.write("KStar: logZLB " + dfp.getLogZ());
                    }

                    fw.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                System.out.println("Epsilon Approx BB: " + logZ);
            }
            if (false) {
                PartFuncTree tree = new PartFuncTree(sp.emat, sp.pruneMat);
                long startTime = System.currentTimeMillis();
                double maxTime = 3600000;
                double logZ = tree.computeEpsilonApprox(0.1, maxTime);
                long totalTime = (System.currentTimeMillis() - startTime);
                System.out.println("New Alg Took: " + totalTime + " milliseconds");
                File statistics = new File("data_var_p_2.txt");
                try {
                    FileWriter fw = new FileWriter(statistics);
                    fw.write("LogConfSpace: " + getLogConfSpace(sp.pruneMat) + "\n");
                    fw.write("NewAlgorithm: " + totalTime + "\n");
                    if (tree.timeOut) {
                        fw.write("NewAlgorithm: finished false" + "\n");
                        fw.write("NewAlgorithm: effectiveEpsilon " + tree.effectiveEpsilon + "\n");
                    } else {
                        fw.write("NewAlgorithm: finished true" + "\n");
                    }
                    fw.write("NewAlgorithm: logZ " + logZ);
                    fw.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
                System.out.println("Epsilon Approx BB: " + logZ);
            }
            /*            ReparamMRF mrf = new ReparamMRF(sp.emat, upm, 0.0);
             MarkovRandomField mrf2 = new MarkovRandomField(sp, 0.0);
             SCMF_Clamp scmf = new SCMF_Clamp(mrf);
             System.out.println("Lower Bound: "+scmf.getLogZLB());
             //            DiscretePartFunc dfp = new DiscretePartFunc(sp.emat, upm, 0.1);
            
             //            TRBP_Refactor_2 trbpR = new TRBP_Refactor_2(mrf2);
             /*            TRBPSeq trbp = new TRBPSeq(mrf);
             double ubLogZ = trbp.getLogZ();
             System.out.println("ubLogZ: " + ubLogZ);
             if (false) {
             DiscretePartFunc dpf = new DiscretePartFunc(sp, 0.1);
             double logZKStar = dpf.getLogZ();
             System.out.println("KStar LogZ: " + logZKStar);
             }
             */
        }
    }

    public VariationalPartFunc(ConfigFileParser aCFP, double epsilon, double eCut, boolean useTRBPSplit) {
        this.cfp = aCFP;
        this.cfp.params.setValue("STERICTHRESH", "1000");
        SearchProblem[] spList = cfp.getMSDSearchProblems();

        for (SearchProblem searchProb : spList) {
            loadEMatandPrune(searchProb, Double.POSITIVE_INFINITY);
        }

        sp = spList[0];

        UpdatedPruningMatrix upm = new UpdatedPruningMatrix(sp.pruneMat);
        prune(sp, upm, 30);
        this.epsilon = epsilon;
        PartFuncTree tree = new PartFuncTree(sp.emat, upm, eCut, useTRBPSplit);
        long startTime = System.currentTimeMillis();
        double maxTime = 3600000;
        double logZ = tree.computeEpsilonApprox(epsilon, maxTime);
        long totalTime = (System.currentTimeMillis() - startTime);
        System.out.println("New Alg Took: " + totalTime + " milliseconds");
        String filename = "data_";
        filename += epsilon + "_1Hour.txt";
        File statistics = new File(filename);
        try {
            FileWriter fw = new FileWriter(statistics);
            fw.write("LogConfSpace: " + getLogConfSpace(upm) + "\n");
            fw.write("NewAlgorithm: " + totalTime + "\n");
            if (tree.timeOut) {
                fw.write("NewAlgorithm: finished false" + "\n");
                fw.write("NewAlgorithm: effectiveEpsilon " + tree.effectiveEpsilon + "\n");
            } else {
                fw.write("NewAlgorithm: finished true" + "\n");
            }
            fw.write("NewAlgorithm: logZ " + logZ + "\n");
            DiscretePartFunc dfp = new DiscretePartFunc(sp.emat, upm, epsilon, maxTime);
            if (dfp.finishedInTime) {
                fw.write("KStar: finished true" + "\n");
                fw.write("KStar: totalTime " + dfp.totalTime);
            } else {
                fw.write("KStar: finished false" + "\n");
                fw.write("KStar: effectiveEpsilon " + dfp.effectiveEpsilonReached + "\n");
                fw.write("KStar: logZLB " + dfp.getLogZ());
            }

            fw.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        System.out.println("Epsilon Approx BB: " + logZ);
    }

    private double getLogConfSpace(PruningMatrix pruneMat) {
        double logConfSpace = 0;
        for (int pos = 0; pos < pruneMat.numPos(); pos++) {
            logConfSpace += Math.log(pruneMat.unprunedRCsAtPos(pos).size());
        }
        return logConfSpace;
    }

    private void prune(SearchProblem sp, PruningMatrix pruneMat, double pruningInterval) {
        Pruner dee = new Pruner(sp, pruneMat, false, Double.POSITIVE_INFINITY, pruningInterval, false, false, false);
        dee.prune("GOLDSTEIN");
    }

    //Loads energy matrices and prune 
    private void loadEMatandPrune(SearchProblem searchProb, double pruningInterval) {
        System.out.println("Precomputing Energy Matrix for " + searchProb.name + " state");
        searchProb.loadEnergyMatrix();

        System.out.println("Initializing Pruning for " + searchProb.name + " state");
        initializePruning(searchProb);
        PruningControl pruning = cfp.setupPruning(searchProb, pruningInterval, false, false);
        pruning.prune();
    }

    private void prune(SearchProblem searchProb, double pruningInterval) {
        PruningControl pruning = cfp.setupPruning(searchProb, pruningInterval, false, false);
        pruning.prune();
    }

    private void initializePruning(SearchProblem searchProblem) {
        //Creates an efficient competitor pruning matrix
        searchProblem.competitorPruneMat = null;
        //prune with 0 interval, anything that survives will be added as a competitor
        PruningControl compPruning = cfp.setupPruning(searchProblem, Double.POSITIVE_INFINITY, false, false);
        compPruning.setOnlyGoldstein(true);
        compPruning.prune();
        searchProblem.competitorPruneMat = searchProblem.pruneMat;
        searchProblem.pruneMat = null;
    }

    private double getLogZGumbelSample(SearchProblem searchProb, int numSamples) {
        double average = 0.0;
        int averageNodesExpanded = 0;
        for (int i = 0; i < numSamples; i++) {
            GumbelMapTree tree = new GumbelMapTree(searchProb);
            tree.nextConf();
            double score = tree.currentBestFeasibleScore;
            average += score;
            double logZestimate = -average / (this.constRT * (i + 1));
            System.out.println("Number of Nodes Expanded: " + tree.numExpanded);
            averageNodesExpanded += tree.numExpanded;
            System.out.println("Current Sample: " + -score / this.constRT);
            System.out.println("Current Average: " + logZestimate);
            System.out.println("Current Average Nodes Exp: " + averageNodesExpanded / (i + 1));
        }
        System.out.println("Average Nodes Expanded: " + averageNodesExpanded / numSamples);
        double logZ = -average / (this.constRT * numSamples);
        System.out.println("Gumbel logZ: " + logZ);

        return logZ;
    }

    private double kstarScoreGumbelSample(SearchProblem[] spList, int numSamples) {
        double logZBound = getLogZGumbelSample(spList[0], numSamples);
        double logZProtein = getLogZGumbelSample(spList[1], numSamples);
        double logZLigand = getLogZGumbelSample(spList[2], numSamples);

        return logZBound - (logZLigand + logZProtein);
    }

    private HashMap<Integer, Integer> getUnboundPosNum2BoundPosNum(SearchProblem unboundSP, SearchProblem boundSP) {
        List<String> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<String> resNumsBound = boundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        HashMap<Integer, Integer> unboundPosToBoundPos = new HashMap<>();
        for (int pos = 0; pos < unboundSP.confSpace.numPos; pos++) {
            //Get the index of the corresponding resNum for the unbound pos
            //in the boun resNum list.
            unboundPosToBoundPos.put(pos, resNumsBound.indexOf(resNumsUnbound.get(pos)));
        }

        return unboundPosToBoundPos;
    }

    private void getBoundPosNum2UnboundPosNum(SearchProblem unboundSP, SearchProblem boundSP, HashMap<Integer, Integer> bound2Unbound) {
        List<String> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<String> resNumsBound = boundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        for (int pos = 0; pos < unboundSP.confSpace.numPos; pos++) {
            //Get the index of the corresponding resNum for the unbound pos
            //in the boun resNum list.
            if (resNumsUnbound.contains(resNumsBound.get(pos))) {
                bound2Unbound.put(pos, resNumsBound.indexOf(resNumsUnbound.get(pos)));
            }
        }
    }

    private void getBoundPosNum2UnboundSearchProblem(SearchProblem unboundSP, SearchProblem boundSP, HashMap<Integer, SearchProblem> bound2Unbound) {
        List<String> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<String> resNumsBound = boundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        for (int pos = 0; pos < unboundSP.confSpace.numPos; pos++) {
            //Get the index of the corresponding resNum for the unbound pos
            //in the boun resNum list.
            if (resNumsUnbound.contains(resNumsBound.get(pos))) {
                bound2Unbound.put(pos, unboundSP);
            }
        }
    }

    private int[][] getConfSamples(SearchProblem searchProb, int numSamples) {
        int[][] samples = new int[numSamples][searchProb.confSpace.numPos];

        for (int i = 0; i < numSamples; i++) {
            GumbelMapTree tree = new GumbelMapTree(searchProb);
            tree.nextConf();
            samples[i] = tree.currentBestFeasibleSolution;
        }

        return samples;
    }

    private int[][] getBoundSamplesFromUnboundSamples(SearchProblem boundSP, SearchProblem unboundSP1, int[][] samples1, SearchProblem unboundSP2, int[][] samples2) {
        //For now I will just combine the samples additively instead of multiplicatively 
        if (samples1.length != samples2.length) {
            throw new RuntimeException("Error: Unbound States Have Different Number of Samples");
        }

        HashMap<Integer, Integer> unboundSP1Pos2BoundPos = getUnboundPosNum2BoundPosNum(unboundSP1, boundSP);
        HashMap<Integer, Integer> unboundSP2Pos2BoundPos = getUnboundPosNum2BoundPosNum(unboundSP2, boundSP);

        int numSamples = samples1.length;
        int[][] boundSamples = new int[numSamples][boundSP.confSpace.numPos];

        for (int i = 0; i < numSamples; i++) {
            int[] unbound1Sample = samples1[i];
            int[] unbound2Sample = samples2[i];

            int[] boundSample = new int[boundSP.confSpace.numPos];
            for (int pos1 = 0; pos1 < unboundSP1.confSpace.numPos; pos1++) {
                int boundPosNum = unboundSP1Pos2BoundPos.get(pos1);
                boundSample[boundPosNum] = unbound1Sample[pos1];
            }
            for (int pos2 = 0; pos2 < unboundSP2.confSpace.numPos; pos2++) {
                int boundPosNum = unboundSP2Pos2BoundPos.get(pos2);
                boundSample[boundPosNum] = unbound2Sample[pos2];
            }

            boundSamples[i] = boundSample;
        }

        return boundSamples;
    }

    private double[] getSampleEnergies(SearchProblem searchProb, int[][] samples) {
        int numSamples = samples.length;

        double[] sampleEs = new double[numSamples];
        for (int i = 0; i < numSamples; i++) {
            int[] conf = samples[i];
            double E = searchProb.emat.getInternalEnergy(new RCTuple(conf));
            sampleEs[i] = E;
        }

        return sampleEs;
    }

    private double getLogKstarSampling(SearchProblem[] spList, int numSamples) {
        SearchProblem boundSP = spList[0];
        SearchProblem unboundSP1 = spList[1];
        SearchProblem unboundSP2 = spList[2];

        int[][] samplesUnbound1 = getConfSamples(unboundSP1, numSamples);
        double[] sampleEnergyUnbound1 = getSampleEnergies(unboundSP1, samplesUnbound1);
        int[][] samplesUnbound2 = getConfSamples(unboundSP2, numSamples);
        double[] sampleEnergyUnbound2 = getSampleEnergies(unboundSP2, samplesUnbound2);
        int[][] samplesBound = getBoundSamplesFromUnboundSamples(boundSP, unboundSP1, samplesUnbound1, unboundSP2, samplesUnbound2);
        double[] sampleEnergyBound = getSampleEnergies(boundSP, samplesBound);

        double averageDiffE = computeAverageDiffEnergy(sampleEnergyBound, sampleEnergyUnbound1, sampleEnergyUnbound2);
        return -averageDiffE / this.constRT;
    }

    private double getLogKstarSampling2(SearchProblem[] spList, int numSamples) {
        SearchProblem boundSP = spList[0];
        SearchProblem unboundSP1 = spList[1];
        SearchProblem unboundSP2 = spList[2];

        int[][] samplesBound = getConfSamples(boundSP, numSamples);
        double[] boundSamplesE = getSampleEnergies(boundSP, samplesBound);
        HashMap<Integer, SearchProblem> bound2UnboundSearchProb = new HashMap<>();
        getBoundPosNum2UnboundSearchProblem(unboundSP1, boundSP, bound2UnboundSearchProb);
        getBoundPosNum2UnboundSearchProblem(unboundSP2, boundSP, bound2UnboundSearchProb);

        int[][] unboundSample1 = getUnboundSamples(unboundSP1, samplesBound, bound2UnboundSearchProb);
        double[] unboundSamples1E = getSampleEnergies(unboundSP1, unboundSample1);
        int[][] unboundSample2 = getUnboundSamples(unboundSP2, unboundSample1, bound2UnboundSearchProb);
        double[] unboundSamples2E = getSampleEnergies(unboundSP2, unboundSample2);

        double averageDiffE = computeAverageDiffEnergy(boundSamplesE, unboundSamples1E, unboundSamples2E);
        return averageDiffE / this.constRT;
    }

    private int[][] getUnboundSamples(SearchProblem unboundSP, int[][] boundSamples, HashMap<Integer, SearchProblem> bound2UnboundSP) {
        int numSamples = boundSamples.length;
        int[][] unboundSamples = new int[numSamples][unboundSP.confSpace.numPos];

        for (int i = 0; i < numSamples; i++) {
            int[] boundSample = boundSamples[i];
            int[] unboundSample = new int[unboundSP.confSpace.numPos];
            int iter = 0;
            for (int pos = 0; pos < boundSample.length; pos++) {
                if (bound2UnboundSP.get(pos) == unboundSP) {
                    unboundSample[iter] = boundSample[pos];
                    iter++;
                }
            }
            unboundSamples[i] = unboundSample;
        }
        return unboundSamples;
    }

    private double computeAverageDiffEnergy(double[] boundSample, double[] unboundSample1, double[] unboundSample2) {
        if ((boundSample.length != unboundSample1.length) || (boundSample.length != unboundSample2.length)) {
            throw new RuntimeException("Bound and Unbound Samples Have Different Length");
        }
        int numSamples = boundSample.length;

        double average = 0;
        for (int i = 0; i < numSamples; i++) {
            double boundE = boundSample[i];
            double unboundE = unboundSample1[i] + unboundSample2[i];
            average += (boundE - unboundE);
        }
        average = average / numSamples;

        return average;
    }

    private void testSyntheticMST() {
        int numNodes = 5;
        boolean[][] interactionGraph = getCompleteGraph(numNodes);
        double[][] edgeWeights = new double[numNodes][];
        for (int i = 0; i < numNodes; i++) {
            edgeWeights[i] = new double[i];
            for (int j = 0; j < i; j++) {
                if (i - j == 1) {
                    edgeWeights[i][j] = -10;
                } else {
                    edgeWeights[i][j] = 10;
                }
            }
        }
        MinSpanningTree mst = new MinSpanningTree(edgeWeights, interactionGraph);
        double[][] vector = mst.mstVector;
    }

    boolean[][] getCompleteGraph(int numNodes) {
        boolean[][] interactionGraph = new boolean[numNodes][];
        for (int i = 0; i < numNodes; i++) {
            interactionGraph[i] = new boolean[numNodes];
            for (int j = 0; j < numNodes; j++) {
                if (j != i) {
                    interactionGraph[i][j] = true;
                } else {
                    interactionGraph[i][j] = false;
                }
            }
        }
        return interactionGraph;
    }

    private double computeGMECLB(SearchProblem sp) {
        UpdatedPruningMatrix pm = new UpdatedPruningMatrix(sp.pruneMat);
        Pruner dee = new Pruner(sp, pm, false, Double.POSITIVE_INFINITY,
                0.0, false, false, false);
        dee.prune("GOLDSTEIN");
        ConfTree tree = new ConfTree(sp.emat, pm);
        int[] conf = tree.nextConf();
        double gmecE = sp.emat.getInternalEnergy(new RCTuple(conf));
        return -gmecE / constRT;
    }

}
