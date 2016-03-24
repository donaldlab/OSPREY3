/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.kadee.GumbelMapTree;
import edu.duke.cs.osprey.astar.partfunc.partFuncTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.partitionfunctionbounds.MapPerturbation;
import edu.duke.cs.osprey.partitionfunctionbounds.MarkovRandomField;
import edu.duke.cs.osprey.partitionfunctionbounds.SelfConsistentMeanField;
import edu.duke.cs.osprey.partitionfunctionbounds.TreeReweightedBeliefPropagation;
import edu.duke.cs.osprey.pruning.PruningControl;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 *
 * @author hmn5
 */
public class GumbelDoer {

    ConfigFileParser cfp;
    SearchProblem sp;
    final double constRT = PoissonBoltzmannEnergy.constRT;

    public GumbelDoer(ConfigFileParser aCFP) {
        this.cfp = aCFP;
        SearchProblem[] spList = cfp.getMSDSearchProblems();

        for (SearchProblem searchProb : spList) {
            loadEMatandPrune(searchProb, Double.POSITIVE_INFINITY);
        }

        SearchProblem sp = spList[0];
        partFuncTree tree = new partFuncTree(sp);
        double Z = tree.computeEpsilonApprox(0.1);
        
        /*
        MarkovRandomField mrf = new MarkovRandomField(sp, 0.0);
        SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
        scmf.run();
        double lb = scmf.calcLBLogZ();
        System.out.println("Lower Bound: "+lb);
        //TreeReweightedBeliefPropagation trbp = new TreeReweightedBeliefPropagation(mrf);
        MapPerturbation mp = new MapPerturbation(sp);
        int numBelow = 0;
        double average = 0.0;
        for (int i = 0; i < 1000; i++) {
            double ub = mp.calcUBLogZLPMax(10);
            if (ub < lb) {
                System.out.println("Upper: " + ub + "  Lower: " + lb);
                numBelow++;
            }
            average += ub;
            System.out.println("MPLP Iter: " + i + "  Sample: " + ub);
        }
        System.out.println("Num Below: " + numBelow);
        System.out.println("Average: "+average/1000);
        double kstarSCMF = kstarScoreSCMF(spList);
        System.out.println("KStar SCMF: " + kstarSCMF);
        System.out.println();

//        double kstarGumbelSample = kstarScoreGumbelSample(spList, 300);
        double kstarSample = getLogKstarSampling2(spList, 300);
        System.out.println("KStar GumbelSample: " + kstarSample);
*/

        /*
         sp = cfp.getSearchProblem();
         loadEMatandPrune(this.sp, Double.POSITIVE_INFINITY);

         double average = 0.0;
         int averageNodesExpanded = 0;
         int numSamples = 1;
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

         MarkovRandomField mrf = new MarkovRandomField(sp, 0.0);
         SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
         scmf.run();
         double lowerZ = scmf.calcLBLogZ();
         System.out.println("SCMF LogZ: " + lowerZ);
         MapPerurbation mpert = new MapPerurbation(sp, true);
         double upperZ = mpert.calcUBLogZ(200);
         System.out.println("MapPert LogZ:" + upperZ); */
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

    private void initializePruning(SearchProblem searchProblem) {
        //Creates an efficient competitor pruning matrix
        searchProblem.competitorPruneMat = null;
        System.out.println("PRECOMPUTING COMPETITOR PRUNING MATRIX");
        //prune with 0 interval, anything that survives will be added as a competitor
        PruningControl compPruning = cfp.setupPruning(searchProblem, Double.POSITIVE_INFINITY, false, false);
        compPruning.setOnlyGoldstein(true);
        compPruning.prune();
        searchProblem.competitorPruneMat = searchProblem.pruneMat;
        searchProblem.pruneMat = null;
        System.out.println("COMPETITOR PRUNING DONE");
    }

    private double kstarScoreSCMF(SearchProblem[] spList) {
        double logZBound = getLogZSCMF(spList[0]);
        System.out.println("logZBound SCMF: " + logZBound);
        double logZProtein = getLogZSCMF(spList[1]);
        System.out.println("logZProtein SCMF: " + logZProtein);
        double logZLigand = getLogZSCMF(spList[2]);
        System.out.println("logZLigand SCMF: " + logZLigand);

        return logZBound - (logZLigand + logZProtein);
    }

    private double getLogZSCMF(SearchProblem searchProb) {
        MarkovRandomField mrf = new MarkovRandomField(searchProb, 0.0);
        SelfConsistentMeanField scmf = new SelfConsistentMeanField(mrf);
        scmf.run();
        double logZ = scmf.calcLBLogZ();

        return logZ;
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
        List<Integer> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<Integer> resNumsBound = boundSP.confSpace.posFlex.stream()
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
        List<Integer> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<Integer> resNumsBound = boundSP.confSpace.posFlex.stream()
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
        List<Integer> resNumsUnbound = unboundSP.confSpace.posFlex.stream()
                .map(posFlex -> posFlex.res.resNum)
                .collect(Collectors.toCollection(ArrayList::new));

        List<Integer> resNumsBound = boundSP.confSpace.posFlex.stream()
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
}
