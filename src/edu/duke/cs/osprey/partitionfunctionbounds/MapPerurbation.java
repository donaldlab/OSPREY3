/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.ObjectIO;
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class MapPerurbation {

    final SearchProblem searchSpace;
    EnergyMatrix emat;

    final double constRT = PoissonBoltzmannEnergy.constRT;
    final ExpFunction ef = new ExpFunction();

    //Keep track of the confs of the MAP that we get for each sample
    ///We can analyze these with mutual information to determine how to merge
    public int[][] mapConfsUB;
    public int[][] mapConfsLB;

    //For analyzing probabilities
    int numSamplesAnalysis;
    ArrayList<singlePos> singlePosList;
    ArrayList<pairPos> pairPosList;

    GumbelDistribution gd; 
    public MapPerurbation(SearchProblem searchSpace) {
        this.searchSpace = searchSpace;
        this.emat = searchSpace.emat;
        gd = new GumbelDistribution();
    }

    //Returns Upper Bounds on Log Partition Function
    public double calcUBLogZ(int anumSamples) {
        int numSamples = anumSamples;
        mapConfsUB = new int[numSamples][emat.oneBody.size()];
        BigDecimal averageGMECs = new BigDecimal(0.0);

        for (int i = 0; i < numSamples; i++) {
            ArrayList<ArrayList<Double>> originalOneBodyEmat = (ArrayList<ArrayList<Double>>) ObjectIO.deepCopy(emat.oneBody);
            addUBGumbelNoiseOneBody();
            ConfSearch search = new ConfTree(searchSpace);
            int[] conf = search.nextConf();
            mapConfsUB[i] = conf;
            double E = -1.0 * searchSpace.lowerBound(conf);
            averageGMECs = averageGMECs.add(new BigDecimal(E));
            //replace oneBody with original to remove the noise added
            emat.oneBody = originalOneBodyEmat;
        }

        return averageGMECs.divide(new BigDecimal(numSamples * this.constRT), ef.mc).doubleValue();
    }

    //Returns Lower Bound on natural log of Partition Function
    public double calcLBLogZ(int anumSamples) {
        int numSamples = anumSamples;
        mapConfsLB = new int[numSamples][emat.oneBody.size()];
        BigDecimal averageGMECs = new BigDecimal(0.0);

        for (int i = 0; i < numSamples; i++) {
            ArrayList<ArrayList<Double>> originalOneBodyEmat = (ArrayList<ArrayList<Double>>) ObjectIO.deepCopy(emat.oneBody);
            addLBGumbelNoiseOneBody();
            ConfSearch search = new ConfTree(searchSpace);
            int[] conf = search.nextConf();
            mapConfsLB[i] = conf;
            double E = -1.0 * searchSpace.lowerBound(conf);
            averageGMECs = averageGMECs.add(new BigDecimal(E));
            //replace oneBody with original to remove the noise added
            emat.oneBody = originalOneBodyEmat;
            if (i%100==0){
                System.out.println("Map Pert Iteration: "+i);
            }
        }

        return averageGMECs.divide(new BigDecimal(numSamples * this.constRT), ef.mc).doubleValue();
    }

    //Returns lower bound on log_10 of partition function
    public double calcLBLog10Z(int aNumSamples) {
        return (Math.log10(Math.E)) * calcLBLogZ(aNumSamples);
    }
    
    public double calcUBLog10Z(int aNumSamples){
        return (Math.log10(Math.E))*calcUBLogZ(aNumSamples);
    }


    //add Gumbel noise to one-body terms
    private void addUBGumbelNoiseOneBody() {
        for (int pos = 0; pos < emat.oneBody.size(); pos++) {
            for (int superRC : searchSpace.pruneMat.unprunedRCsAtPos(pos)) {
                double currentE = emat.getOneBody(pos, superRC);
                double noise = gd.sample(-1.0 * GumbelDistribution.gamma, 1.0) * this.constRT;
                emat.setOneBody(pos, superRC, currentE - noise);
            }
        }
    }

    //add Gumbel noise to one-body terms
    private void addLBGumbelNoiseOneBody() {
        for (int pos = 0; pos < emat.oneBody.size(); pos++) {
            for (int superRC : searchSpace.pruneMat.unprunedRCsAtPos(pos)) {
                double currentE = emat.getOneBody(pos, superRC);
                double noise = gd.sample(-1.0 * GumbelDistribution.gamma, 1.0) * this.constRT / emat.oneBody.size();
                emat.setOneBody(pos, superRC, currentE - noise);
            }
        }
    }

    //get the counts of each rotamer and each pair of rotamers at each pos or pair of pos
    public void getMapPosRotCounts(int[][] mapConfs) {
        this.numSamplesAnalysis = mapConfs.length;
        int numPos = mapConfs[0].length;

        this.singlePosList = new ArrayList<>();
        this.pairPosList = new ArrayList<>();

        for (int sample = 0; sample < numSamplesAnalysis; sample++) {
            int[] conf = mapConfs[sample];
            for (int posNum1 = 0; posNum1 < numPos; posNum1++) {
                //Get res
                if (singlePosList.size() <= posNum1) {
                    singlePos res = new singlePos(posNum1);
                    singlePosList.add(res);
                }
                singlePos pos1 = singlePosList.get(posNum1);
                //Get rot
                int rotNum1 = conf[posNum1];
                singleRot rot1 = new singleRot(rotNum1);
                //contains uses equals which is overriden 
                if (!pos1.rotList.contains(rot1)) {
                    pos1.rotList.add(rot1);
                } else {
                    //update count
                    int indexRot1 = pos1.rotList.indexOf(rot1);
                    pos1.rotList.get(indexRot1).count++;
                }

                for (int posNum2 = 0; posNum2 < posNum1; posNum2++) {
                    int rotNum2 = conf[posNum2];
                    pairPos twoPos = new pairPos(posNum1, posNum2);
                    if (!pairPosList.contains(twoPos)) {
                        pairPosList.add(twoPos);
                    }
                    int indexPairPos = pairPosList.indexOf(twoPos);
                    pairPos pair = pairPosList.get(indexPairPos);
                    pairRot rotPair = new pairRot(rotNum1, rotNum2);
                    if (!pair.rotList.contains(rotPair)) {
                        pair.rotList.add(rotPair);
                    } else {
                        int indexRotPair = pair.rotList.indexOf(rotPair);
                        pair.rotList.get(indexRotPair).count++;
                    }
                }
            }
        }
    }

    public ArrayList<Integer> getPairWithMaxMutualInfo(boolean useUpperBoundMapConfs) {
        //returns arraylist of size 2, where each element is a position in the pair
        //with maximum mutual information
        if (useUpperBoundMapConfs) {
            getMapPosRotCounts(mapConfsUB);
        } else {
            getMapPosRotCounts(mapConfsLB);
        }
        double maxMutualInfo = 0.0;
        pairPos maxPair = new pairPos(-1, -1);
        for (pairPos pair : pairPosList) {
            double mutualInfo = calcMutualInformation(pair);
            if (mutualInfo > maxMutualInfo) {
                maxMutualInfo = mutualInfo;
                maxPair = pair;
            }
        }
        ArrayList<Integer> pairPosNums = new ArrayList<>();
        if (maxPair.res1 < maxPair.res2) {
            pairPosNums.add(maxPair.res1);
            pairPosNums.add(maxPair.res2);
        } else {
            pairPosNums.add(maxPair.res2);
            pairPosNums.add(maxPair.res1);
        }
        return pairPosNums;
    }

    private double calcMutualInformation(pairPos pair) {
        singlePos pos1 = singlePosList.get(pair.res1);
        singlePos pos2 = singlePosList.get(pair.res2);
        double mutualInfo = calcSinglePosEntropy(pos1) + calcSinglePosEntropy(pos2) - calcPairPosEntropy(pair);
        return mutualInfo;
    }

    private double calcSinglePosEntropy(singlePos pos) {
        double[] probabilities = new double[pos.rotList.size()];
        //get probabilities
        for (int i = 0; i < pos.rotList.size(); i++) {
            singleRot rot = pos.rotList.get(i);
            double prob = (double) rot.count / (double) this.numSamplesAnalysis;
            probabilities[i] = prob;
        }
        //entropy 
        double entropy = 0.0;
        for (double prob : probabilities) {
            entropy += -1.0 * (prob * Math.log(prob));
        }
        return entropy;
    }

    private double calcPairPosEntropy(pairPos pair) {
        double[] probabilities = new double[pair.rotList.size()];
        //get probabilities
        for (int i = 0; i < pair.rotList.size(); i++) {
            pairRot rot = pair.rotList.get(i);
            double prob = (double) rot.count / (double) this.numSamplesAnalysis;
            probabilities[i] = prob;
        }
        //entropy 
        double entropy = 0.0;
        for (double prob : probabilities) {
            entropy += -1.0 * (prob * Math.log(prob));
        }
        return entropy;
    }

    private class singlePos {

        int resNum;
        ArrayList<singleRot> rotList = new ArrayList<>();

        public singlePos(int resNum) {
            this.resNum = resNum;
        }

        @Override
        public boolean equals(Object ares2) {
            singlePos res2 = (singlePos) ares2;
            return this.resNum == res2.resNum;
        }
    }

    private class singleRot {

        int rotNum;
        int count = 1;

        public singleRot(int rotNum) {
            this.rotNum = rotNum;
        }

        @Override
        public boolean equals(Object asingleRot2) {
            singleRot singleRot2 = (singleRot) asingleRot2;
            return this.rotNum == singleRot2.rotNum;
        }
    }

    private class pairPos {

        int res1;
        int res2;

        ArrayList<pairRot> rotList = new ArrayList<>();

        public pairPos(int res1, int res2) {
            this.res1 = res1;
            this.res2 = res2;
        }

        @Override
        public boolean equals(Object aresPair2) {
            pairPos resPair2 = (pairPos) aresPair2;
            boolean match1 = (this.res1 == resPair2.res1) && (this.res2 == resPair2.res2);
            boolean match2 = (this.res1 == resPair2.res2) && (this.res2 == resPair2.res1);
            boolean same = match1 || match2;
            return same;
        }
    }

    private class pairRot {

        int rot1;
        int rot2;
        int count = 1;

        public pairRot(int rot1, int rot2) {
            this.rot1 = rot1;
            this.rot2 = rot2;
        }

        @Override
        public boolean equals(Object apairRot2) {
            pairRot pairRot2 = (pairRot) apairRot2;
            boolean match1 = (this.rot1 == pairRot2.rot1) && (this.rot2 == pairRot2.rot2);
            boolean match2 = (this.rot1 == pairRot2.rot2) && (this.rot2 == pairRot2.rot1);
            return match1 || match2;
        }
    }
}

