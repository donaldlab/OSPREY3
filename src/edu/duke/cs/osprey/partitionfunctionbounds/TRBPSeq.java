/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.CreateMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class TRBPSeq {

    double logZ = Double.POSITIVE_INFINITY;

    //Node list
    ArrayList<MRFNode> nodeList;
    //Emat defines the potential functions
    int numNodes;
    int[] numLabelsPerNode;

    UpdatedEmat emat;
    //interactionGraph;
    boolean[][] nonClampledInteractionGraph;

    TupleMatrix<Double> marginalProbabilities;

    //threshold for convergence
    double threshold = 1e-6;
    final int maxIterations = 50;

    double constRT = PoissonBoltzmannEnergy.constRT;

    double damping = 0.5;
    int numEdgeProbUpdates = 15;

    //p_e are the edge probabilities over spanning trees
    public double[][] edgeProbabilities;
    double[][][] messages;

    double[][] expNormMessages;
    double[][] expNormMarginals;
    boolean dampLog = false;

    public ExpFunction ef = new ExpFunction();

    double maxChange;
    double averageChange;

    boolean verbose = true;
    boolean printDuringEdgeUpdate = true;
    boolean useArmijosRule = false;

    public TRBPSeq(ReparamMRF mrf) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.nonClampledInteractionGraph = mrf.nonClampedInteractionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.labelList.size();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.nonClampledInteractionGraph);
        this.messages = initializeMessages(1.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);

        runTRBPSeq5();
    }

    private void runTRBPSeq() {
        for (int j = 0; j < numEdgeProbUpdates; j++) {
            computeExpNormals();
            int numIter = 0;
            while (maxChange > threshold || numIter < 2) {
                maxChange = 0;
                numIter++;
                double[][][] messagesNPlus1 = updateMessages(this.messages);
//                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
                if (true) {
                    updateMarginals(this.messages);
//                    checkMarginals();
                    double currentlogZ = calcUBLogZ();
                    if (verbose) {
                        System.out.println("Max Change in Messages: " + this.maxChange + "  LogZUB: " + currentlogZ);
                    }
                }
            }
//            System.out.println("TRBP took: " + numIter + " iterations");
            updateMarginals(this.messages);
//            checkMarginals();
            //System.out.println("LogZUB: "+logZ);

            double currentlogZ = calcUBLogZ();
            if (!Double.isNaN(currentlogZ) && !Double.isInfinite(currentlogZ)) {
                this.logZ = Math.min(this.logZ, currentlogZ);
            }
            if (j < numEdgeProbUpdates - 1) {
                MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), nonClampledInteractionGraph);
                double stepSize = getStepSizeParameter(j + 1);
                updateEdgeProbabilies(mst.mstVector, stepSize);
                checkEdgeProbabilities();
            }
        }
    }

    private void runTRBPSeq2() {
        for (int j = 0; j < numEdgeProbUpdates; j++) {
            computeExpNormals();
            int numIter = 0;
            while (((maxChange > threshold) && (numIter < this.maxIterations)) || numIter < 2) {
                maxChange = 0;
                numIter++;
//                double[][][] messagesNPlus1 = updateMessages(this.messages);
                double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
                updateMarginals(this.messages);
//                    checkMarginals();
                double currentlogZ = calcUBLogZ();
                this.logZ = currentlogZ;
                if (verbose) {
                    System.out.println("Max Change in Messages: " + this.maxChange + "  LogZUB: " + currentlogZ);
                }
            }
//            System.out.println("TRBP took: " + numIter + " iterations");
//            updateMarginals(this.messages);
//            checkMarginals();
//            System.out.println("LogZUB: "+logZ);

//            double currentlogZ = calcUBLogZ();
//            System.out.println("LogZUB: " + currentlogZ);
//            if (!Double.isNaN(currentlogZ)) {
//                this.logZ = Math.min(this.logZ, currentlogZ);
//            }
            if (j < numEdgeProbUpdates - 1) {
                System.out.println("logZ UB: " + this.logZ);
                if (verbose) {
                    System.out.println("Updating Edge Probabilities:  ");
                }
                MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), nonClampledInteractionGraph);
                double stepSize = getStepSizeParameter(j + 1);
                updateEdgeProbabilies(mst.mstVector, stepSize);
                checkEdgeProbabilities();
            } else {
                double change = Double.POSITIVE_INFINITY;
                while (change < 1e-2) {
                    double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
                    checkMessages(messagesNPlus1);
                    this.messages = messagesNPlus1;
                    updateMarginals(this.messages);
//                    checkMarginals();
                    double currentlogZ = calcUBLogZ();
                    change = Math.abs(this.logZ - currentlogZ);
                    this.logZ = currentlogZ;
                }
                System.out.println("LogZ UB: " + this.logZ);
            }
        }
    }

    private void runTRBPSeq3() {
        double changeBetweenEdgeUpdates = Double.POSITIVE_INFINITY;
        double lastLogZBetweenUpdates = Double.POSITIVE_INFINITY;
        int numEdgeUpdates = 0;
        double[][] edgeProbGradient = this.edgeProbabilities;
        while (changeBetweenEdgeUpdates > 0.1) {
            double changeWithinEdgeUpdate = Double.POSITIVE_INFINITY;
            double lastLogZ = Double.POSITIVE_INFINITY;
            int numUpdatesWithinEdgeProb = 0;
            //The first update does not change the initial edge probabilites
            //That is why we initizliae the edgeProbGradient with the current edge probabilies.
            //Otherwise, this gradient is a 0-1 vector
            double stepSize = getStepSizeParameter(numEdgeUpdates + 1);
            updateEdgeProbabilies(edgeProbGradient, stepSize);

            checkEdgeProbabilities();
            computeExpNormals();
            this.messages = initializeMessages(1.0);
            while (changeWithinEdgeUpdate > 0.01) {
                double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
//                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
                updateMarginals(this.messages);
                double currentlogZ = calcUBLogZ();
                if (verbose) {
                    System.out.println("    logZ UB: " + currentlogZ);
                }
                changeWithinEdgeUpdate = Math.abs(lastLogZ - currentlogZ);
                lastLogZ = currentlogZ;
                numUpdatesWithinEdgeProb++;
            }
            System.out.println("Updating Edge Probabilites After " + numEdgeUpdates + " iterations. LogZ UB: " + lastLogZ);
            this.logZ = Math.min(this.logZ, lastLogZ);
            changeBetweenEdgeUpdates = Math.abs(lastLogZBetweenUpdates - lastLogZ);
            lastLogZBetweenUpdates = lastLogZ;
            MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), nonClampledInteractionGraph);
            edgeProbGradient = mst.mstVector;
            numEdgeUpdates++;
        }
        System.out.println("Updated Edge Probabilities " + (numEdgeUpdates - 1) + " times. LogZ UB: " + this.logZ);
    }

    private void runTRBPSeq4() {
        double changeBetweenEdgeUpdates = Double.POSITIVE_INFINITY;
        double lastLogZBetweenUpdates = Double.POSITIVE_INFINITY;
        int numEdgeUpdates = 0;
        double[][] edgeProbGradient = this.edgeProbabilities;
        while (numEdgeUpdates < 5) {
            double changeWithinEdgeUpdate = Double.POSITIVE_INFINITY;
            double lastLogZ = Double.POSITIVE_INFINITY;
            int numUpdatesWithinEdgeProb = 0;
            //The first update does not change the initial edge probabilites
            //That is why we initizliae the edgeProbGradient with the current edge probabilies.
            //Otherwise, this gradient is a 0-1 vector
            double stepSize = getStepSizeParameter(numEdgeUpdates + 1);
            updateEdgeProbabilies(edgeProbGradient, stepSize);

            checkEdgeProbabilities();
            computeExpNormals();
            this.messages = initializeMessages(1.0);
            while (numUpdatesWithinEdgeProb < 30) {
                double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
//                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
                updateMarginals(this.messages);
                double currentlogZ = calcUBLogZ();
                if (verbose) {
                    System.out.println("    logZ UB: " + currentlogZ);
                }
                changeWithinEdgeUpdate = Math.abs(lastLogZ - currentlogZ);
                lastLogZ = currentlogZ;
                numUpdatesWithinEdgeProb++;
            }
            System.out.println("Updating Edge Probabilites After " + numEdgeUpdates + " iterations. LogZ UB: " + lastLogZ);
            this.logZ = Math.min(this.logZ, lastLogZ);
            changeBetweenEdgeUpdates = Math.abs(lastLogZBetweenUpdates - lastLogZ);
            lastLogZBetweenUpdates = lastLogZ;
            MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), nonClampledInteractionGraph);
            edgeProbGradient = mst.mstVector;
            numEdgeUpdates++;
        }

        double stepSize = getStepSizeParameter(numEdgeUpdates + 1);
        updateEdgeProbabilies(edgeProbGradient, stepSize);

        checkEdgeProbabilities();
        computeExpNormals();
        this.messages = initializeMessages(1.0);
        double changeWithinEdgeUpdate = Double.POSITIVE_INFINITY;
        double lastLogZ = Double.POSITIVE_INFINITY;
        int numUpdatesWithinEdgeProb = 0;
        while (changeWithinEdgeUpdate > .01) {
            double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
//                checkMessages(messagesNPlus1);
            this.messages = messagesNPlus1;
            updateMarginals(this.messages);
            double currentlogZ = calcUBLogZ();
            if (verbose) {
                System.out.println("    logZ UB: " + currentlogZ);
            }
            changeWithinEdgeUpdate = Math.abs(lastLogZ - currentlogZ);
            lastLogZ = currentlogZ;
            numUpdatesWithinEdgeProb++;
        }
        System.out.println("TRBPS Finished after " + numUpdatesWithinEdgeProb + " iterations   logZUB: " + lastLogZ);
        this.logZ = lastLogZ;
    }

    private void runTRBPSeq5() {
        double changeBetweenEdgeUpdates = Double.POSITIVE_INFINITY;
        double lastLogZBetweenUpdates = Double.POSITIVE_INFINITY;
        int numEdgeUpdates = 0;
        double[][] edgeProbGradient = this.edgeProbabilities;

//        double[][] edgeProb = GraphUtils.getEdgeProbabilities(nonClampledInteractionGraph);
//        this.edgeProbabilities = GraphUtils.getEdgeProbabilities(this.nonClampledInteractionGraph);
        double stepSize = 1.0;
//        computeExpNormals();
        while (changeBetweenEdgeUpdates > 0.01) {
            double changeWithinEdgeUpdate = Double.POSITIVE_INFINITY;
            double lastLogZ = Double.POSITIVE_INFINITY;
            int numUpdatesWithinEdgeProb = 0;
            //The first update does not change the initial edge probabilites
            //That is why we initizliae the edgeProbGradient with the current edge probabilies.
            //Otherwise, this gradient is a 0-1 vector

            this.messages = initializeMessages(1.0);
            if (numEdgeUpdates > 0) {
                updateEdgeProbabilies(edgeProbGradient, stepSize);
            }
            computeExpNormals();
//            checkEdgeProbabilities();
            while (changeWithinEdgeUpdate > 0.001 || numUpdatesWithinEdgeProb < 50) {
                double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
//                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
                updateMarginals(this.messages);
                double currentlogZ = calcUBLogZ();
                if (verbose) {
                    System.out.println("    logZ UB: " + currentlogZ + "   max change: " + this.maxChange + "  average: " + this.averageChange);
                }
                changeWithinEdgeUpdate = Math.abs(lastLogZ - currentlogZ);
                lastLogZ = currentlogZ;
                numUpdatesWithinEdgeProb++;
            }
            System.out.println("Updating Edge Probabilites After " + numEdgeUpdates + " iterations. LogZ UB: " + lastLogZ);
            this.logZ = Math.min(this.logZ, lastLogZ);
            changeBetweenEdgeUpdates = Math.abs(lastLogZBetweenUpdates - lastLogZ);
            lastLogZBetweenUpdates = lastLogZ;
            MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), nonClampledInteractionGraph);
            edgeProbGradient = mst.mstVector;

            stepSize = getStepSizeParameter(numEdgeUpdates + 1);
            if (changeBetweenEdgeUpdates < 0.01) {
                changeWithinEdgeUpdate = Double.POSITIVE_INFINITY;
                while (changeWithinEdgeUpdate > 0.0001) {
                    double[][][] messagesNPlus1 = updateMessagesSeq(this.messages);
//                checkMessages(messagesNPlus1);
                    this.messages = messagesNPlus1;
                    updateMarginals(this.messages);
                    double currentlogZ = calcUBLogZ();
                    if (verbose) {
                        System.out.println("    logZ UB: " + currentlogZ);
                    }
                    changeWithinEdgeUpdate = Math.abs(lastLogZ - currentlogZ);
                    lastLogZ = currentlogZ;
                    numUpdatesWithinEdgeProb++;
                }
                lastLogZBetweenUpdates = lastLogZ;
            }

            numEdgeUpdates++;
        }
        this.logZ = Math.min(this.logZ, lastLogZBetweenUpdates);
        System.out.println("Updated Edge Probabilities " + (numEdgeUpdates - 1) + " times. LogZ UB: " + this.logZ);
    }

    private void computeExpNormals() {
        double[][] expNormMessage = new double[this.numNodes][this.numNodes];
        double[][] expNormMarginal = new double[this.numNodes][this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        MRFLabel rotI = nodeI.labelList.get(rotIindex);
                        MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                        double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                        double edgeProbability = getEdgeProbability(i, j);
                        double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
                        double rotIE = -this.emat.getOneBody(nodeI.nodeNum, rotI.labelNum);
                        double toBeExponentiateJI = (pairwiseE / edgeProbability) + rotJE;
                        double toBeExponentiateIJ = (pairwiseE / edgeProbability) + rotIE;
                        double toBeExponentiatedMarginal = (pairwiseE / edgeProbability) + rotIE + rotJE;
                        expNormMessage[j][i] = Math.max(expNormMessage[j][i], toBeExponentiateJI);
                        expNormMessage[i][j] = Math.max(expNormMessage[i][j], toBeExponentiateIJ);
                        expNormMarginal[i][j] = Math.max(expNormMarginal[i][j], toBeExponentiatedMarginal);
                        expNormMarginal[j][i] = expNormMarginal[i][j];
                    }
                }
            }
        }
        this.expNormMessages = expNormMessage;
        this.expNormMarginals = expNormMarginal;
    }

    private double[][][] updateMessages(double[][][] previousMessages) {
        double[][][] updatedMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);
        double maxChangeMessage = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {

                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                if (this.nonClampledInteractionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        double messageAtRotI = 0.0;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
//                            double normalized = ((pairwiseE / edgeProbability) + rotJE); //- expNormI;
                            double normalized = ((pairwiseE / edgeProbability) + rotJE - this.expNormMessages[j][i]);
                            double messageFromRotJ = Math.exp(normalized / this.constRT);
                            //double messageFromRotJ = Math.exp(normalized);
                            messageFromRotJ = messageFromRotJ * (getProductMessages(nodeJ, rotJindex, nodeI, previousMessages));
                            messageAtRotI = messageAtRotI + (messageFromRotJ);
                        }
                        updatedMessages[j][i][rotIindex] = messageAtRotI;
                        partFunctionI += messageAtRotI;
                    }
                    if (!dampLog) {
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = damping * updatedMessages[j][i][rotIindex] + ((1 - damping) * previousMessages[j][i][rotIindex]);
                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    } else {
                        double partFuncI2 = 0.0;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = Math.pow(updatedMessages[j][i][rotIindex], damping) + Math.pow(previousMessages[j][i][rotIindex], (1 - damping));
                            partFuncI2 += updatedMessages[j][i][rotIindex];
                        }
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFuncI2;
                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    }
                    //Now we send message from i->j
                    double partFunctionJ = 0.0;
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        double messageAtRotJ = 0.0;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);

                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotIE = -this.emat.getOneBody(nodeI.nodeNum, rotI.labelNum);
//                            double normalized = ((pairwiseE / edgeProbability) + rotIE); //- expNormJ;
                            double normalized = ((pairwiseE / edgeProbability) + rotIE - this.expNormMessages[i][j]);
                            double messageFromRotI = Math.exp(normalized / this.constRT);
                            //double messageFromRotI = Math.exp(normalized);
                            messageFromRotI = messageFromRotI * (getProductMessages(nodeI, rotIindex, nodeJ, previousMessages));
                            messageAtRotJ += messageFromRotI;
                        }
                        updatedMessages[i][j][rotJindex] = messageAtRotJ;
                        partFunctionJ += messageAtRotJ;
                    }
                    if (!dampLog) {
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            updatedMessages[i][j][rotJindex] /= partFunctionJ;
                            updatedMessages[i][j][rotJindex] = damping * updatedMessages[i][j][rotJindex] + ((1 - damping) * previousMessages[i][j][rotJindex]);
                            double prevMessage = previousMessages[i][j][rotJindex];
                            double change = Math.abs(prevMessage - updatedMessages[i][j][rotJindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go

                            checkNumericalStability(updatedMessages, i, j, rotJindex);
                            previousMessages[i][j][rotJindex] = updatedMessages[i][j][rotJindex];
                        }
                    } else {
                        double partFuncJ2 = 0.0;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            updatedMessages[i][j][rotJindex] /= partFunctionJ;
                            updatedMessages[i][j][rotJindex] = Math.pow(updatedMessages[i][j][rotJindex], damping) + Math.pow(previousMessages[i][j][rotJindex], (1 - damping));
                            partFuncJ2 += updatedMessages[i][j][rotJindex];
                        }
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            updatedMessages[i][j][rotJindex] /= partFuncJ2;
                            double prevMessage = previousMessages[i][j][rotJindex];
                            double change = Math.abs(prevMessage - updatedMessages[i][j][rotJindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, i, j, rotJindex);
                            previousMessages[i][j][rotJindex] = updatedMessages[i][j][rotJindex];
                        }
                    }
                }
            }
        }
        this.maxChange = maxChangeMessage;
        return updatedMessages;
    }

    private double[][][] updateMessagesSeq(double[][][] previousMessages) {
        double[][][] updatedMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);
        double maxChangeMessage = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {

                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                if (this.nonClampledInteractionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        double messageAtRotI = 0.0;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
//                            double normalized = ((pairwiseE / edgeProbability) + rotJE); //- expNormI;
                            double normalized = ((pairwiseE / edgeProbability) + rotJE - this.expNormMessages[j][i]);
                            double messageFromRotJ = Math.exp(normalized / this.constRT);
                            //double messageFromRotJ = Math.exp(normalized);
                            messageFromRotJ = messageFromRotJ * (getProductMessages(nodeJ, rotJindex, nodeI, previousMessages));
                            messageAtRotI = messageAtRotI + (messageFromRotJ);
                        }
                        updatedMessages[j][i][rotIindex] = messageAtRotI;
                        partFunctionI += messageAtRotI;
                    }
                    if (!dampLog) {
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = damping * updatedMessages[j][i][rotIindex] + ((1 - damping) * previousMessages[j][i][rotIindex]);
                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    } else {
                        double partFuncI2 = 0.0;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = Math.pow(updatedMessages[j][i][rotIindex], damping) + Math.pow(previousMessages[j][i][rotIindex], (1 - damping));
                            partFuncI2 += updatedMessages[j][i][rotIindex];
                        }
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFuncI2;

                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    }
                }
            }
        }
        double averChange = 0.0;
        int numMessages = 0;
        for (int i = this.numNodes - 1; i > -1; i--) {
            for (int j = this.numNodes - 1; j > i; j--) {

                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                if (this.nonClampledInteractionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        double messageAtRotI = 0.0;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
//                            double normalized = ((pairwiseE / edgeProbability) + rotJE); //- expNormI;
                            double normalized = ((pairwiseE / edgeProbability) + rotJE - this.expNormMessages[j][i]);
                            double messageFromRotJ = Math.exp(normalized / this.constRT);
                            //double messageFromRotJ = Math.exp(normalized);
                            messageFromRotJ = messageFromRotJ * (getProductMessages(nodeJ, rotJindex, nodeI, previousMessages));
                            messageAtRotI = messageAtRotI + (messageFromRotJ);
                        }
                        updatedMessages[j][i][rotIindex] = messageAtRotI;
                        partFunctionI += messageAtRotI;
                    }
                    if (!dampLog) {
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = damping * updatedMessages[j][i][rotIindex] + ((1.0 - damping) * previousMessages[j][i][rotIindex]);
                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                            averChange += change;
                            numMessages++;
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    } else {
                        double partFuncI2 = 0.0;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFunctionI;
                            updatedMessages[j][i][rotIindex] = Math.pow(updatedMessages[j][i][rotIindex], damping) + Math.pow(previousMessages[j][i][rotIindex], (1 - damping));
                            partFuncI2 += updatedMessages[j][i][rotIindex];
                        }
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            updatedMessages[j][i][rotIindex] /= partFuncI2;
                            double prevMessage = previousMessages[j][i][rotIindex];
                            double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);

                            averChange += change;
                            numMessages++;
                            maxChangeMessage = Math.max(change, maxChangeMessage);
                            //Update as we go
                            checkNumericalStability(updatedMessages, j, i, rotIindex);
                            previousMessages[j][i][rotIindex] = updatedMessages[j][i][rotIindex];
                        }
                    }
                }
            }
        }

        this.maxChange = maxChangeMessage;
        this.averageChange = averChange / numMessages;
        return updatedMessages;
    }

    private void checkNumericalStability(double[][][] aMessage, int sendingNode, int receivingNode, int rot) {
        if (aMessage[sendingNode][receivingNode][rot] < 1e-100) {
            aMessage[sendingNode][receivingNode][rot] = 1e-100;
        }
    }

    private void checkMessages(double[][][] messages) {
        double minMessage = Double.POSITIVE_INFINITY;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                double sumJ = 0.0;
                for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                    if (messages[i][j][rotJ] == 0.0) {
                        throw new RuntimeException("MESSAGE = 0");
                    }
                    sumJ += messages[i][j][rotJ];
                    minMessage = Math.min(minMessage, messages[i][j][rotJ]);
                }
                if (Math.abs(1 - sumJ) > 1e-6) {
                    throw new RuntimeException("Messages Not Normalized between node " + nodeI.nodeNum + " and node " + nodeJ.nodeNum);
                }

                double sumI = 0.0;
                for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                    if (messages[j][i][rotI] == 0.0) {
                        throw new RuntimeException("MESSAGE = 0");
                    }
                    sumI += messages[j][i][rotI];
                    minMessage = Math.min(minMessage, messages[j][i][rotI]);
                }
                if (Math.abs(1 - sumI) > 1e-6) {
                    throw new RuntimeException("Messages Not Normalized between node " + nodeJ.nodeNum + " and node " + nodeI.nodeNum);
                }
            }
        }
    }

    private void updateMarginals(double[][][] messages) {
        //first update node marginals
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            updateNodeMarginal(node, messages);
        }
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (this.nonClampledInteractionGraph[i][j]) {
                    MRFNode nodeI = this.nodeList.get(i);
                    MRFNode nodeJ = this.nodeList.get(j);
                    updateEdgeMarginal(nodeI, nodeJ, messages);
                }
            }
        }}

    

    private void updateEdgeMarginal(MRFNode nodeI, MRFNode nodeJ, double[][][] messages) {
        double partFuncOverall = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double pairwiseE = this.emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double edgeProb = getEdgeProbability(nodeI.index, nodeJ.index);
                double nodeIE = this.emat.getOneBody(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum);
                double nodeJE = this.emat.getOneBody(nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);

//                double normalized = (-((pairwiseE / edgeProb) + nodeIE + nodeJE)); //- expNorm;
                double normalized = (-((pairwiseE / edgeProb) + nodeIE + nodeJE)) - this.expNormMarginals[nodeI.index][nodeJ.index];
                double marginal = Math.exp(normalized / this.constRT);
//                double marginal = Math.exp(normalized);
                marginal = marginal * (getProductMessages(nodeI, rotI, nodeJ, messages));
                marginal = marginal * (getProductMessages(nodeJ, rotJ, nodeI, messages));

                partFuncOverall += marginal;
                this.marginalProbabilities.setPairwise(nodeI.index, rotI, nodeJ.index, rotJ, marginal);
            }
        }
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double unNormalized = this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
                double normalized = unNormalized / partFuncOverall;
                this.marginalProbabilities.setPairwise(nodeI.index, rotI, nodeJ.index, rotJ, normalized);
            }
        }
    }

    private void checkMarginals() {
        for (int i = 0; i < this.numNodes; i++) {
            checkNodeMarginal(this.nodeList.get(i));
            for (int j = 0; j < i; j++) {
                checkPairwiseMarginal(this.nodeList.get(i), this.nodeList.get(j));
            }
        }
    }

    private void checkNodeMarginal(MRFNode node) {
        double sum = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            /*            if (this.marginalProbabilies.getOneBody(node.index, rot) == 0.0){
             throw new RuntimeException("Margine = 0");
             }*/
            sum += this.marginalProbabilities.getOneBody(node.index, rot);
        }
        if (Math.abs(sum - 1) > 1e-6) {
            System.out.println("Error in Marginal: " + Math.abs(sum - 1));
            System.out.println("Marginal Not Normalized at Node with PosNum: " + node.nodeNum);
        }
    }

    private void checkPairwiseMarginal(MRFNode nodeI, MRFNode nodeJ) {
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            double marginal = this.marginalProbabilities.getOneBody(nodeI.index, rotI);
            double sum = 0.0;
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                sum += this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
            }
            if (Math.abs(sum - marginal) > 0.0001) {
                System.out.println("Error in Marginal: " + Math.abs(sum - marginal));
                System.out.println("Marginal Not Normalized at Node Pairs: " + nodeI.nodeNum + " " + nodeJ.nodeNum);
            }
        }
        for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
            double marginal = this.marginalProbabilities.getOneBody(nodeJ.index, rotJ);
            double sum = 0.0;
            for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                sum += this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
            }
            if (Math.abs(sum - marginal) > 0.0001) {
                System.out.println("Error in Marginal: " + Math.abs(sum - marginal));
                System.out.println("Marginal Not Normalized at Node Pairs with Pos Nums: " + nodeJ.nodeNum + " " + nodeI.nodeNum);
            }
        }
    }

    private void updateNodeMarginal(MRFNode node, double[][][] messages) {
        double partFunc = 0.0;
        double expNorm = Double.NEGATIVE_INFINITY;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            expNorm = Math.max(expNorm, -this.emat.getOneBody(node.nodeNum, node.labelList.get(rot).labelNum));
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double normalized = -this.emat.getOneBody(node.nodeNum, node.labelList.get(rot).labelNum) - expNorm;
            double update = Math.exp(normalized / this.constRT);
//            double update = Math.exp(normalized);
            for (MRFNode neighbor : node.neighborList) {
                double message = messages[neighbor.index][node.index][rot];
                double prob = getEdgeProbability(neighbor.index, node.index);
                update = update * (Math.pow(message, prob));
            }
            this.marginalProbabilities.setOneBody(node.index, rot, update);
            partFunc += update;
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double unNormalized = this.marginalProbabilities.getOneBody(node.index, rot);
            this.marginalProbabilities.setOneBody(node.index, rot, unNormalized / partFunc);
        }
    }

    private double getProductMessages(MRFNode nodeReceiving, int labelIndex, MRFNode nodeExcluded, double[][][] currentMessages) {
        double product = 1;
        for (MRFNode node : nodeReceiving.neighborList) {
            if (!(node.equals(nodeExcluded))) {
                double message = currentMessages[node.index][nodeReceiving.index][labelIndex];
                double probability = getEdgeProbability(nodeReceiving.index, node.index);
                product = product * (Math.pow(message, probability));
            }
        }
        double messageNodeExcluded = currentMessages[nodeExcluded.index][nodeReceiving.index][labelIndex];
        double probability = getEdgeProbability(nodeExcluded.index, nodeReceiving.index);
        product = product / (Math.pow(messageNodeExcluded, (1 - probability)));
        return product;
    }

    private double getEdgeProbability(int nodeI, int nodeJ) {
        if (nodeI < nodeJ) {
            return this.edgeProbabilities[nodeJ][nodeI];
        }
        return this.edgeProbabilities[nodeI][nodeJ];
    }

    private double[][][] initializeMessages(double initVal) {
        if (initVal <= 0) {
            throw new RuntimeException("TRBP ERROR: Initial Message Value Must be Positive");
        }

        double[][][] messages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, initVal);
        //TODO NORMALIZE!!!!!!!!!!!
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                //First normalize i-j x_j
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                    messages[nodeI.index][nodeJ.index][rotJ] = 1.0 / (double) nodeJ.labelList.size();
                }
                //Then normalize j-i x_i
                for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                    messages[nodeJ.index][nodeI.index][rotI] = 1.0 / (double) nodeI.labelList.size();
                }
            }
        }
        return messages;
    }

    /**
     * We initialize the edge probabilities by 1/numEdges for each edge in the
     * graph
     *
     * @param interactionGraph
     * @return
     */
    private double[][] initializeEdgeProbabilities(boolean[][] interactionGraph) {
        int numEdges = 0;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (interactionGraph[i][j]) {
                    numEdges++;
                }
            }
        }

        double[][] edgeProb = new double[this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            //to avoid overcounting, we index between neighbors i,j where j<i
            double[] edgesFromI = new double[i];
            for (int j = 0; j < i; j++) {
                if (interactionGraph[i][j]) {
                    //TODO: This is a uniform distribution over a complete graph
                    //This should be generalized to non-complete graphs
                    edgesFromI[j] = 2.0 / ((double) this.nodeList.size());
                } else {
                    edgesFromI[j] = 0.0;
                }
            }
            edgeProb[i] = edgesFromI;
        }
        return edgeProb;
    }

    private double getSingleNodeEnthalpy(MRFNode node) {
        double enthalpy = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double E = this.emat.getOneBody(node.nodeNum, node.labelList.get(rot).labelNum);
            double prob = this.marginalProbabilities.getOneBody(node.index, rot);
            enthalpy += E * prob;
        }
        return enthalpy;
    }

    private double getPairwiseNodeEnthalpy(MRFNode nodeI, MRFNode nodeJ) {
        double enthalpy = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double E = emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double prob = this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
                enthalpy += E * prob;
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
                if (this.nonClampledInteractionGraph[node1.index][node2.index]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double prob = this.marginalProbabilities.getOneBody(node.index, rot);
            if (prob != 0.0) {
                double entropyAtRot = (-1.0) * prob * Math.log(prob);
                if (Double.isFinite(entropyAtRot)) {
                    entropy += entropyAtRot;
                }
            }
        }
        return entropy;
    }

    private double getMutualInformation(MRFNode nodeI, MRFNode nodeJ) {
        double mutualInf = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double probIJ = this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
                double probI = this.marginalProbabilities.getOneBody(nodeI.index, rotI);
                double probJ = this.marginalProbabilities.getOneBody(nodeJ.index, rotJ);
//                if (probIJ != 0.0) {
                if ((probIJ != 0.0) && (probI != 0.0) && (probJ != 0.0)) {
                    double mutualInfAtRotPair = probIJ * Math.log(probIJ / (probI * probJ));
                    if (Double.isFinite(mutualInfAtRotPair)) {
                        mutualInf += probIJ * Math.log(probIJ / (probI * probJ));
                    }
                }
            }
        }
        return mutualInf;
    }

    /**
     * Returns true if any probabilities are zero
     *
     * @param probabilities
     * @return
     */
    private boolean probabilitiesEqualZero(double... probabilities) {
        for (double prob : probabilities) {
            if (prob == 0.0) {
                return true;
            }
        }
        return false;
    }

    public double getEntropy() {
        double entropy = 0.0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode nodeI = this.nodeList.get(i);
            entropy += getSingleNodeEntropy(nodeI);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (nonClampledInteractionGraph[i][j]) {
                    double edgeProb = getEdgeProbability(i, j);
                    entropy -= edgeProb * getMutualInformation(nodeI, nodeJ);
                }
            }
        }
        return entropy;
    }

    private double calcFreeEnergy() {
        double enthalpy = getEnthalpy();
        double entropy = getEntropy();
        double freeEnergy = enthalpy - this.constRT * entropy;
        return freeEnergy;
    }

    public double calcUBLogZ() {
        return -(calcFreeEnergy() + this.emat.getConstTerm()) / this.constRT;
    }

    double[][] getEdgeWeights() {
        double[][] edgeWeights = new double[this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            double[] edgeFromI = new double[i];
            MRFNode nodeI = this.nodeList.get(i);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (nonClampledInteractionGraph[i][j]) {
                    edgeFromI[j] = -getMutualInformation(nodeI, nodeJ);
                }
            }
            edgeWeights[i] = edgeFromI;
        }
        return edgeWeights;
    }

    private double getStepSizeParameter(int iteration) {
        return 2.0 / (iteration + 4.0);
    }

    private void updateEdgeProbabilies(double[][] descentDirection, double stepSize) {
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                this.edgeProbabilities[i][j] = stepSize * descentDirection[i][j] + (1 - stepSize) * edgeProbabilities[i][j];
            }
        }
    }

    public static double round(double value, int places) {
        if (places < 0) {
            throw new IllegalArgumentException();
        }

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

    private double[][] copy2DArray(double[][] original) {
        double[][] copy = new double[original.length][];
        for (int i = 0; i < original.length; i++) {
            copy[i] = original[i].clone();
        }
        return copy;
    }

    public double getLogZ() {
        return this.logZ;
    }

    private void checkEdgeProbabilities() {
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                double prob = this.edgeProbabilities[i][j];
                if (prob > 1.0) {
                    throw new RuntimeException("Edge Probability is Greater Than 1: " + prob);
                }
            }
        }
    }
}
