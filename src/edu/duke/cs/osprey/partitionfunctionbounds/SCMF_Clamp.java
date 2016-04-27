/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class SCMF_Clamp {

    //Nodelist defines the graphical model, since nodes keep track of their neighbors
    ArrayList<MRFNode> nodeList;
    int numNodes;
    //Emat defines the potential functions for the Probabilistic Graphical Model
    UpdatedEmat emat;
    //interactionGraph determines if two nodes are neighbors
    boolean[][] nonClampedInteractionGraph;

    //to avoid numerical overflow we computer "exp normalizers"
    double[] expNormMessages;

    //threshold for convergence (change in beliefs must be less than threshold)
    double threshold = 1e-8;
    final int minNumberIterations = 300;
    final int maxNumberIterations = 10000;
    double lambda = 1.0;

    double scmfTemp = PoissonBoltzmannEnergy.constRT;
    public ExpFunction ef = new ExpFunction();

    double freeEnergy;
    boolean hasCalculatedFreeEnergy;

    boolean verbose = false;
    boolean debugMode = false;

    public SCMF_Clamp(ReparamMRF mrf) {
        this.nodeList = mrf.nodeList;
        this.numNodes = mrf.nodeList.size();
        this.emat = mrf.emat;
        this.nonClampedInteractionGraph = mrf.nonClampedInteractionGraph;
//        double startTime = System.currentTimeMillis();
        run();
//        double totalTime = System.currentTimeMillis() - startTime;
//        System.out.println("Took "+totalTime+" milliseconds to run");
    }

    public void run() {
        int iter = 0;
        boolean hasConverged = false;
        initializeBeliefs();
        //to avoid numerical overflow we computer "exp normalizers"
        computeExpNormals();
        //set the temperature to be high
//        this.scmfTemp = this.scmfTemp;
        //while we haven't converged and are below maxNumberIteration;
        while (!hasConverged && (iter < this.maxNumberIterations)) {
            //update all beliefs

            hasConverged = updateAllBeliefs();

            if (iter == 100) {
                this.freeEnergy = calcFreeEnergy();
                this.scmfTemp = this.scmfTemp * 1000;
//                lambda = 1.0;
            }
            if (debugMode) {
                checkBeliefs();
            }
            //lower the temperature
            lowerTemperature(0.98);
            //Make sure we are at least beyond the minimum number of iterations
            hasConverged = hasConverged && (iter > this.minNumberIterations);
            iter++;
        }
        this.freeEnergy = Math.min(this.freeEnergy, calcFreeEnergy());
        if (verbose) {
            printMaxLabel();
        }
        resetTemperature();
    }
    
    private void run2() {
        int iter = 0;
        int numIterBeforeTempRaise = 500;
        boolean hasConverged = false;
        initializeBeliefs();
        //to avoid numerical overflow we computer "exp normalizers"
        computeExpNormals();
        //set temperature to be high
        this.scmfTemp = this.scmfTemp * 1000;
        while (!hasConverged && (iter < maxNumberIterations)) {
            //update all beliefs
            hasConverged = updateAllBeliefs();
            if (debugMode) {
                checkBeliefs();
            }
            //lower the temperature
            double tempMultiplier = 0.98;
            lowerTemperature(tempMultiplier);
            if (iter == numIterBeforeTempRaise) {
                //store this free energy in cases it is better than after raising temperature
                this.freeEnergy = calcFreeEnergy();
                hasCalculatedFreeEnergy = true;
                this.scmfTemp = this.scmfTemp * 1000;
            } else if (iter == numIterBeforeTempRaise * 2) {
                //take best free energy so far
                this.freeEnergy = Math.min(this.freeEnergy, calcFreeEnergy());
                hasCalculatedFreeEnergy = true;
                this.scmfTemp = this.scmfTemp * 500;
            }
            //Make sure we are at least beyond the minimum number of iterations
            hasConverged = hasConverged && (iter > this.minNumberIterations);
            iter++;
        }
        //take best free energy
        this.freeEnergy = Math.min(this.freeEnergy, calcFreeEnergy());
        hasCalculatedFreeEnergy = true;
        if (verbose) {
            System.out.println("SCMF Finished After " + iter + " iterations");
            printMaxLabel();
        }
        resetTemperature();
    }

    //updates all beliefs and returns true if 
    //max difference in beliefs is less than threshold
    private boolean updateAllBeliefs() {
        boolean allNodesConverged = true;
        for (MRFNode node : this.nodeList) {
            //We converge if and only if every node has converged
            boolean nodeConverged = updateNodeBeliefs(node);
            if (!(allNodesConverged && nodeConverged)) {
                allNodesConverged = false;
            }
        }
        return allNodesConverged;
    }

    //update all of the beliefs belonging to node
    //returns true if max difference in beleifs is less than threshold
    private boolean updateNodeBeliefs(MRFNode node) {
        //create a normalizing constant to normalize beliefs
        double partFunction = 0.0;
        //keep track of difference between beliefs for convergence;
        double maxEpsilon = 0.0;

        int numLabels = node.labelList.size();
        double[] updatedBeliefs = new double[numLabels];

        //iterate over labels to get partition function value
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            double oneBodyE = this.emat.getOneBody(node.posNum, label.labelNum);
            double meanFieldE = getMeanFieldEnergy(node, label);
            //unnormalized updateBelief
            double logUnnormalizedBelief = -(oneBodyE + meanFieldE) / scmfTemp;
            double logExpNormalizedBelief = logUnnormalizedBelief - (this.expNormMessages[node.index] / scmfTemp);
            updatedBeliefs[labelIndex] = Math.exp(logExpNormalizedBelief);
            partFunction += Math.exp(logExpNormalizedBelief);
        }
        //now we update the beliefs using our partFunction (normalizing constant)
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            //keep track of old belief
            double oldBelief = label.currentBelief;
            //if partFunction is 0.0, all labels had really bad energies and we will make them uniform
            if (partFunction == 0.0) {
                label.currentBelief = 1.0 / node.labelList.size();
            } //otherwise we update the beliefs normally
            else {
                label.currentBelief = updatedBeliefs[labelIndex] / partFunction;
                label.currentBelief = lambda * label.currentBelief + (1.0 - lambda) * oldBelief;
            }
            //update maxEpsilon
            if (Math.abs(label.currentBelief - oldBelief) > maxEpsilon) {
                maxEpsilon = Math.abs(label.currentBelief - oldBelief);
            }
        }
        boolean hasConverged = false;
        if (maxEpsilon < threshold) {
            hasConverged = true;
        }
        return hasConverged;
    }

    //lower the temperature, never getting below the true value;
    private void lowerTemperature(double multiplier) {
        this.scmfTemp = this.scmfTemp * multiplier;
        if (this.scmfTemp < PoissonBoltzmannEnergy.constRT) {
            this.scmfTemp = PoissonBoltzmannEnergy.constRT;
        }
    }

    private void resetTemperature() {
        this.scmfTemp = PoissonBoltzmannEnergy.constRT;
    }

    private double getSingleNodeEnthalpy(MRFNode node) {
        double enthalpy = 0.0;
        for (MRFLabel label : node.labelList) {
            double E = this.emat.getOneBody(node.posNum, label.labelNum);
            enthalpy += E * label.currentBelief;
        }
        return enthalpy;
    }

    private double getPairwiseNodeEnthalpy(MRFNode node1, MRFNode node2) {
        double enthalpy = 0.0;
        for (MRFLabel label1 : node1.labelList) {
            for (MRFLabel label2 : node2.labelList) {
                double E = emat.getPairwise(node1.posNum, label1.labelNum, node2.posNum, label2.labelNum);
                enthalpy += E * label1.currentBelief * label2.currentBelief;
            }
        }
        return enthalpy;
    }

    private double getEnthalpy() {
        double enthalpy = 0.0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode node1 = nodeList.get(i);
            enthalpy += getSingleNodeEnthalpy(node1);
            for (int j = 0; j < i; j++) {
                MRFNode node2 = nodeList.get(j);
                if (this.nonClampedInteractionGraph[node1.index][node2.index]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (MRFLabel label : node.labelList) {
            if (label.currentBelief == 0.0) {
                entropy += 0.0;
            } else {
                entropy += (-1.0) * label.currentBelief * Math.log(label.currentBelief);
            }
        }
        return entropy;
    }

    private double getEntropy() {
        double entropy = 0.0;
        for (MRFNode node : this.nodeList) {
            entropy += getSingleNodeEntropy(node);
        }
        return entropy;
    }

    private double calcFreeEnergy() {
        double enthalpy = getEnthalpy();
        double entropy = getEntropy();

        double freeE = enthalpy - this.scmfTemp * entropy;

        return freeE;
    }

    public double getLogZLB() {
        return -(this.freeEnergy + this.emat.getConstTerm()) / this.scmfTemp;
    }

    //Returns the MeanFieldEnergy for a label
    private double getMeanFieldEnergy(MRFNode node, MRFLabel label) {
        double meanFieldE = 0.0;
        for (MRFNode neighbor : node.neighborList) {
            //get the average energy between neighbor and our label
            double averageE = getAverageEnergy(node, label, neighbor);
            meanFieldE += averageE;
        }
        return meanFieldE;
    }

    //Get average energy of label1 from node1, with respect to all labels of node2
    //this is used as a subroutine in getMeanField()
    private double getAverageEnergy(MRFNode node1, MRFLabel label1, MRFNode node2) {
        double averageE = 0.0;
        for (MRFLabel label2 : node2.labelList) {
            double E = this.emat.getPairwise(node1.posNum, label1.labelNum, node2.posNum, label2.labelNum);
            averageE += E * label2.currentBelief;
        }
        return averageE;
    }

    private void initializeBeliefs() {
        for (MRFNode node : this.nodeList) {
            int numLabels = node.labelList.size();
            for (MRFLabel label : node.labelList) {
                label.currentBelief = 1.0 / (double) numLabels;
                //for parallel updates
                label.oldBelief = 1.0 / (double) numLabels;
            }
        }
    }

    public boolean checkBeliefs() {
        boolean beliefsGood = true;
        for (MRFNode node : this.nodeList) {
            double beliefSum = 0.0;
            for (MRFLabel label : node.labelList) {
                beliefSum += label.currentBelief;
            }
            if (!equals(beliefSum, 1.0, 0.001)) {
                beliefsGood = false;
            }
        }
        return beliefsGood;
    }

    public void printMaxLabel() {
        String maxLabels = "";
        for (MRFNode node : nodeList) {
            double max = 0.0;
            int maxLabelNum = -1;
            for (MRFLabel label : node.labelList) {
                if (label.currentBelief > max) {
                    max = label.currentBelief;
                    maxLabelNum = label.labelNum;
                }
            }
            maxLabels = maxLabels + maxLabelNum + " ";
        }
        System.out.println(maxLabels);
    }

    public boolean equals(double a, double b, double epsilon) {
        return a == b ? true : Math.abs(a - b) < epsilon;
    }

    private void computeExpNormals() {
        this.expNormMessages = new double[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            double minOneBody = getMinOneBodyE(node);
            double minMeanField = getMinMeanFieldE(node);
            expNormMessages[i] = -(minOneBody + minMeanField);
        }
    }

    /**
     * returns the minimum possible mean field energy for a node this is the sum
     * of the min over all possible pairwise energies with the node's neighbor
     *
     * @param node
     * @return
     */
    private double getMinMeanFieldE(MRFNode node) {
        double minMeanFieldE = 0.0;
        for (MRFNode neighbor : node.neighborList) {
            double minPairwiseE = Double.POSITIVE_INFINITY;
            for (MRFLabel nodeLabel : node.labelList) {
                for (MRFLabel neighborLabel : neighbor.labelList) {
                    double pairwiseE = this.emat.getPairwise(node.posNum, nodeLabel.labelNum, neighbor.posNum, neighborLabel.labelNum);
                    minPairwiseE = Math.min(minPairwiseE, pairwiseE);
                }
            }
            minMeanFieldE += minPairwiseE;
        }
        return minMeanFieldE;
    }

    private double getMinOneBodyE(MRFNode node) {
        double minE = Double.POSITIVE_INFINITY;
        for (MRFLabel label : node.labelList) {
            double oneBodyE = this.emat.getOneBody(node.posNum, label.labelNum);
            minE = Math.min(minE, oneBodyE);
        }
        return minE;
    }

}
