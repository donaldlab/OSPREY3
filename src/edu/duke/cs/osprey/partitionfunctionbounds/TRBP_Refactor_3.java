/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.confspace.TupleMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.PoissonBoltzmannEnergy;
import edu.duke.cs.osprey.tools.CreateMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class TRBP_Refactor_3 {

    double logZ = Double.POSITIVE_INFINITY;

    ArrayList<MRFNode> nodeList;
    int numNodes;
    int[] numLabelsPerNode;

    EnergyMatrix emat;

    boolean[][] interactionGraph;

    TupleMatrix<Double> marginalProbabilities;

    double threshold = 1e-6;
    double constRT = PoissonBoltzmannEnergy.constRT;

    double damping = 0.5;
    static int maxNumEdgeUpdates = 3;

    public double[][] edgeProbabilities;
    double[][] edgeWeights;
    double[][] edgeClampWeights;

    double[][][] logMessages;
//    double[][][] messages;
    int numMessages;

    boolean useLogDomain = true;
    boolean debug = true;

    public ExpFunction ef = new ExpFunction();

    double maxChange;
    double averageChange;

    double accuracyWithinEdgeProb = 0.001;
    double accuracyBetweenEdgeProb = 0.001;

    double parentUpperBound = Double.POSITIVE_INFINITY;
    boolean useParentUpperBound = false;
    double cutOffUpperBound = 3;

    boolean verbose = false;

    public double[] nodeWeights;

    public TRBP_Refactor_3(MarkovRandomField mrf) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.interactionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.labelList.size();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.interactionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);
        this.numMessages = 2 * getNumEdges(this.interactionGraph);

        initializeEdgeWeights();

        runTRBP();
    }

    public TRBP_Refactor_3(MarkovRandomField mrf, double parentUpperBound) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.interactionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.labelList.size();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.interactionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);
        this.numMessages = 2 * getNumEdges(this.interactionGraph);

        initializeEdgeWeights();

        this.parentUpperBound = parentUpperBound;
        this.useParentUpperBound = true;

        runTRBP();
    }

    private void runTRBP() {
        int numEdgeUpdates = 0;

        //Keep track of last logZ from previous edge update
        double lastLogZEdge = Double.POSITIVE_INFINITY;

        while (numEdgeUpdates <= maxNumEdgeUpdates) {
            if (numEdgeUpdates > 0) {
                if (verbose) {
                    System.out.println("Updating Edge Probabilities...  logZ: " + lastLogZEdge);
                }
                updateEdgeProbabilities(numEdgeUpdates);
                //this.logMessages = initializeLogMessages(0.0);
            }

            double changeBetweenMessageUpdates = Double.POSITIVE_INFINITY;
            //Keep track of last logZ from previous message update
            double lastLogZMessage = Double.POSITIVE_INFINITY;
            int numMessageUpdates = 0;
            //Increase accuracy after initial edge update
            if (numEdgeUpdates == 1) {
                accuracyWithinEdgeProb /= 10.;
            }
            // While we have not converged...
            while ((changeBetweenMessageUpdates > accuracyWithinEdgeProb) || numMessageUpdates < 10) {
                //This is the meat of the algorithm
                updateMessagesParallel(logMessages, damping);
//                updateMessagesSequentially(logMessages, useDamping);
                // Update the marginals and compute LogZ every 10 runs
                if (numMessageUpdates % 10 == 0) {
                    updateMarginals();
                    double currentlogZ = calcUBLogZ();
                    // Get the change in logZ
                    // Classically, people use change in messages rather than change
                    // in logZ, which may be better...
                    double change = Math.abs(lastLogZMessage - currentlogZ);
                    lastLogZMessage = currentlogZ;
                    changeBetweenMessageUpdates = change;
                    if (verbose) {
                        System.out.println("   LogZUB: " + currentlogZ + "  Change: " + change + "  AveMess: " + this.averageChange);
                    }
                }
                numMessageUpdates++;
            }
            lastLogZEdge = lastLogZMessage;
            this.logZ = Math.min(this.logZ, lastLogZEdge);
            // This is (only) useful for the partFuncTree calculation
            // Basically, we don't need high accuracy if the partition function is very small
            // Compared to the parent node 
            if (lastLogZMessage + cutOffUpperBound < parentUpperBound) {
                break;
            }
            numEdgeUpdates++;
        }
        double[] degrees = GraphUtils.getWeightedDegrees(edgeClampWeights);
        this.nodeWeights = degrees;
    }

    int getNumEdges(boolean[][] interactionGraph) {
        int numEdges = 0;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (interactionGraph[i][j]) {
                    numEdges++;
                }
            }
        }
        return numEdges;
    }

    double[][][] initializeLogMessages(double initVal
    ) {
        return CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, initVal);
    }

    double[][] initializeEdgeProbabilities(boolean[][] interactionGraph
    ) {
        return GraphUtils.getEdgeProbabilities(interactionGraph);
    }

    private double getEdgeProbability(MRFNode nodeI, MRFNode nodeJ) {
        int indexNodeI = nodeI.index;
        int indexNodeJ = nodeJ.index;
        if (indexNodeI < indexNodeJ) {
            return this.edgeProbabilities[indexNodeJ][indexNodeI];
        }
        return this.edgeProbabilities[indexNodeI][indexNodeJ];
    }

    private void updateMessagesParallel(double[][][] previousMessages, double damping) {

        double[][][] updatedLogMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);

        int[][] messageOrdering = getNaiveMessagePassingOrdering();
//        int[][] messageOrdering = getMessagePassingOrdering(false);
        int previousSendingNodeIndex = -1;
        double[] sumLogMessages = new double[0];
        for (int[] nodePair : messageOrdering) {
            //Node pair consists of two ints. The first indexes to the sending node
            //The second indexes to the receiving node

            MRFNode nodeI = getNode(nodePair[0]);
            MRFNode nodeJ = getNode(nodePair[1]);

            // Optimization: if we use a particular ordering we can precompute
            // the sumLogMessages for sendingNode and reuse this, which saves a 
            // lot of time
            // If the sending node has changed, we compute the sumLogMessages
            if (previousSendingNodeIndex != nodeI.index) {
                sumLogMessages = getSumLogMessagesV2(nodeI);
                previousSendingNodeIndex = nodeI.index;
            }

            //Get updated log Messages from nodeI to nodeJ
            double[] logMessagesNodeIToNodeJ = getUpdatedLogMessage(nodeI, nodeJ, sumLogMessages);

            double[] previousLogMessages = previousMessages[nodeI.index][nodeJ.index];

            double change = dampenMessages(logMessagesNodeIToNodeJ, previousLogMessages, nodeI.index, nodeJ.index, updatedLogMessages, damping);
            //          this.averageChange += change;
        }

        //Update average change
//        this.averageChange = averageChangeInMessage / this.numMessages;
        this.logMessages = updatedLogMessages;
    }

    private MRFNode getNode(int index) {
        return this.nodeList.get(index);
    }

    public double dampenMessages(double[] updated, double[] previous, int sendIndex, int receiveIndex, double[][][] updatedLogMessages, double damper) {
        double sumChange = 0;
        for (int i = 0; i < updated.length; i++) {
            double prevMessage = previous[i];
            double damped = updated[i] * damper + ((1 - damper) * previous[i]);
            updatedLogMessages[sendIndex][receiveIndex][i] = damped;
            double change = Math.abs(damped - prevMessage);
            sumChange += change;
        }
        return sumChange;
    }

    private void updateMarginals() {
        //first update node marginals
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            updateNodeMarginal(node);
        }
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (this.interactionGraph[i][j]) {
                    MRFNode nodeI = this.nodeList.get(i);
                    MRFNode nodeJ = this.nodeList.get(j);
                    updateEdgeMarginal(nodeI, nodeJ);
                }
            }
        }
    }

    private void updateNodeMarginal(MRFNode node) {
        double partFunc = 0;
        double normalizer = Double.NEGATIVE_INFINITY;
        double[] marginals = new double[node.labelList.size()];
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);

            double labelPot = getOneBodyPotential(node, label) / this.constRT;
//            double normalizedPot = (labelPot - this.expNormNodeMarginals[node.index]) / this.constRT;

            double sumWeightedLogMessages = 0.0;
            for (MRFNode neighbor : this.nodeList) {
                if (this.interactionGraph[node.index][neighbor.index]) {
                    double edgeProb = getEdgeProbability(node, neighbor);
                    double logMessageNeighborToNode = getLogMessage(neighbor, node, label);
                    sumWeightedLogMessages += edgeProb * logMessageNeighborToNode;
                }
            }
            normalizer = Math.max(normalizer, labelPot + sumWeightedLogMessages);
            marginals[labelIndex] = labelPot + sumWeightedLogMessages;
        }
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            marginals[labelIndex] -= normalizer;
            marginals[labelIndex] = Math.exp(marginals[labelIndex]);
            partFunc += marginals[labelIndex];
        }
        for (int labelIndex = 0; labelIndex < node.labelList.size(); labelIndex++) {
            MRFLabel label = node.labelList.get(labelIndex);
            double nonNormalizedMarginal = marginals[labelIndex];
            this.marginalProbabilities.setOneBody(node.index, labelIndex, nonNormalizedMarginal / partFunc);
        }
    }

    private void updateEdgeMarginal(MRFNode nodeI, MRFNode nodeJ) {
        double partFunc = 0;
        double normalizer = Double.NEGATIVE_INFINITY;
        double[] sumLogMessagesItoJ = getSumLogMessages(nodeI, nodeJ);
        double[] sumLogMessagesJtoI = getSumLogMessages(nodeJ, nodeI);
        for (int labelIndexI = 0; labelIndexI < nodeI.labelList.size(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.labelList.size(); labelIndexJ++) {
                MRFLabel labelI = nodeI.labelList.get(labelIndexI);
                MRFLabel labelJ = nodeJ.labelList.get(labelIndexJ);
                double pairPot = getPairwisePotential(nodeI, labelI, nodeJ, labelJ);
                double edgeProb = getEdgeProbability(nodeI, nodeJ);
                double labelIPot = getOneBodyPotential(nodeI, labelI);
                double labelJPot = getOneBodyPotential(nodeJ, labelJ);

                double potential = ((pairPot / edgeProb) + labelIPot + labelJPot) / this.constRT;
//                double normalizedPot = ((pairPot / edgeProb) + labelIPot + labelJPot - this.expNormPairMarginals[nodeI.index][nodeJ.index]) / this.constRT;
//                double sumWeightedLogMessages = getSumLogMessage(nodeI, labelI, nodeJ) + getSumLogMessage(nodeJ, labelJ, nodeI);
                double sumWeightedLogMessages = sumLogMessagesItoJ[labelIndexI + 1] + sumLogMessagesJtoI[labelIndexJ + 1];
                normalizer = Math.max(normalizer, potential + sumWeightedLogMessages);

//                double nonNormalizedLogMarginal = normalizedPot + sumWeightedLogMessages;
//                double nonNormalizedMarginal = Math.exp(nonNormalizedLogMarginal);
//                partFunc += nonNormalizedMarginal;
                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, potential + sumWeightedLogMessages);
            }
        }
        for (int labelIndexI = 0; labelIndexI < nodeI.labelList.size(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.labelList.size(); labelIndexJ++) {
                double unNormalized = this.marginalProbabilities.getPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ);
                double normalized = unNormalized - normalizer;
                double exponentiated = Math.exp(normalized);
                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, exponentiated);
                partFunc += exponentiated;
            }
        }
        for (int labelIndexI = 0; labelIndexI < nodeI.labelList.size(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.labelList.size(); labelIndexJ++) {
                double unNormalized = this.marginalProbabilities.getPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ);
                double normalized = unNormalized / partFunc;
                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, normalized);
            }
        }
    }

    /**
     * Returns a message ordering
     *
     * @return a message ordering of the form
     * int[messageNum][senderNum,receiverNum]
     */
    int[][] getMessagePassingOrdering(boolean sequential) {
        int[][] messagePassingOrdering = new int[numMessages][];

        int currentMessage = 0;
        if (sequential) {
            //First send messages by going from node 1 to node N and sending all 
            //messages backwards (ie. node 2 sends messages to nodes 1 and 0)
            for (int i = 1; i < this.numNodes; i++) {
                for (int j = i - 1; j >= 0; j--) {
                    if (interactionGraph[i][j]) {
                        messagePassingOrdering[currentMessage] = new int[2];
                        messagePassingOrdering[currentMessage][0] = i;
                        messagePassingOrdering[currentMessage][1] = j;
                        currentMessage++;
                    }
                }
            }
            //Now we send messages by going from nodes N-1 to 0 and sending all
            //messages forwards (i.e. node N-2 sends messages to nodes N-1 and N)
            for (int i = this.numNodes - 2; i >= 0; i--) {
                for (int j = i + 1; j < this.numNodes; j++) {
                    if (interactionGraph[i][j]) {
                        messagePassingOrdering[currentMessage] = new int[2];
                        messagePassingOrdering[currentMessage][0] = i;
                        messagePassingOrdering[currentMessage][1] = j;
                        currentMessage++;
                    }
                }
            }
        } else {
            for (int i = 0; i < this.numNodes; i++) {
                for (int j = 0; j < i; j++) {
                    if (interactionGraph[i][j]) {
                        messagePassingOrdering[currentMessage] = new int[2];
                        messagePassingOrdering[currentMessage][0] = i;
                        messagePassingOrdering[currentMessage][1] = j;
                        currentMessage++;
                        messagePassingOrdering[currentMessage] = new int[2];
                        messagePassingOrdering[currentMessage][0] = j;
                        messagePassingOrdering[currentMessage][1] = i;
                        currentMessage++;
                    }
                }
            }
        }
        return messagePassingOrdering;
    }

    //Just for comparison/testing purposes
    int[][] getNaiveMessagePassingOrdering() {
        int[][] messagePassingOrdering = new int[numMessages][];

        int currentMessage = 0;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < this.numNodes; j++) {
                if (interactionGraph[i][j]) {
                    messagePassingOrdering[currentMessage] = new int[2];
                    messagePassingOrdering[currentMessage][0] = i;
                    messagePassingOrdering[currentMessage][1] = j;
                    currentMessage++;
                }
            }
        }
        return messagePassingOrdering;
    }

    double[] getUpdatedLogMessage(MRFNode sendingNode, MRFNode receivingNode, double[] sumLogMessages) {
        double[] updatedLogMessages = new double[receivingNode.getNumLabels()];
        double largestLogMessage = Double.NEGATIVE_INFINITY;
        
        ArrayList<Double> oneBodyESending = this.emat.getOneBodyEnergies(sendingNode.posNum);
        ArrayList<ArrayList<Double>> twoBodyE = this.emat.getTwoBodyEnergies(sendingNode.posNum, receivingNode.posNum);
        
        /*double[] sumLogMessagesMinusReceiver = new double[sendingNode.getNumLabels()];
         for (int i = 0; i < sumLogMessages.length; i++) {
         sumLogMessagesMinusReceiver[i] = sumLogMessages[i] - getLogMessage(receivingNode, sendingNode, sendingNode.labelList.get(i));
         }*/
        double[] sumLogMessagesMinusReceiver = subtractReceiver(sumLogMessages, sendingNode, receivingNode);

        for (MRFLabel receivingLabel : receivingNode.labelList) {
            double normalizer = Double.NEGATIVE_INFINITY;
            double[] toBeSummed = new double[sendingNode.getNumLabels()];

            for (int sendingLabelIndex = 0; sendingLabelIndex < sendingNode.getNumLabels(); sendingLabelIndex++) {
                MRFLabel sendingLabel = sendingNode.labelList.get(sendingLabelIndex);
//                RCTuple tup = new RCTuple(sendingNode.posNum, sendingLabel.labelNum,
//                        receivingNode.posNum, receivingLabel.labelNum);
//                if (pruneMat.isPruned(tup)) {
//                    toBeSummed[sendingLabelIndex] = Double.NEGATIVE_INFINITY;
//                    continue;
//                }
                //compute everything between the curly brackets in eq 39. of TRBP Paper
                //We do this in the log domain and exponentiate aftern
//                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel);
                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel, twoBodyE);
                double edgeProbability = getEdgeProbability(receivingNode, sendingNode);
//                double sendingLabelPot = getOneBodyPotential(sendingNode, sendingLabel);
                double sendingLabelPot = getOneBodyPotential(sendingLabel, oneBodyESending);
                
//                double normalized = ((pairwisePot / edgeProbability) + sendingLabelPot - getExpNormMessagesAtRot(sendingNode, receivingNode, index));
//                double sumLogMessage = sumLogMessages[sendingLabelIndex] - getLogMessage(receivingNode, sendingNode, sendingLabel);
                double sumLogMessage = sumLogMessagesMinusReceiver[sendingLabelIndex];
//                toBeSummed[sendingLabelIndex] = (((pairwisePot / edgeProbability) + sendingLabelPot) / constRT) + sumLogMessage;
                toBeSummed[sendingLabelIndex] = getMessageUpdateInSum(pairwisePot, edgeProbability, sendingLabelPot, sumLogMessage);
                normalizer = Math.max(normalizer, toBeSummed[sendingLabelIndex]);
            }
//            double newNormalizer = expNormMessage/constRT + maxSumLogMessage;
//            normalizer = expNormMessage/constRT + maxSumLogMessage;
            /*double sum = 0.;
             for (int sendingLabelIndex = 0; sendingLabelIndex < sendingNode.labelList.size(); sendingLabelIndex++) {
             toBeSummed[sendingLabelIndex] -= normalizer;
             sum += Math.exp(toBeSummed[sendingLabelIndex]);
             }*/
            int numToBeSummed = sendingNode.getNumLabels();
            double sum = sumNormalizedToBeSummed(toBeSummed, normalizer, numToBeSummed);

//            double nonNormalizedLogMessage = Math.log(sum) + normalizer;
            double nonNormalizedLogMessage = getNonNormalizedLogMessage(sum, normalizer);

            largestLogMessage = Math.max(largestLogMessage, nonNormalizedLogMessage);
            int messageIndex = getLabelIndex(receivingNode, receivingLabel);
            //int messageIndex = receivingNode.labelList.indexOf(receivingLabel);
            updatedLogMessages[messageIndex] = nonNormalizedLogMessage;
        }
        //Normalize by subtracting the largest logMessage;
/*        for (int i = 0; i < updatedLogMessages.length; i++) {
            updatedLogMessages[i] -= largestLogMessage;
        }*/
        normalizeLogMessages(updatedLogMessages, largestLogMessage);
        /*        double partFunc = 0.0;
         for (int i = 0; i < updatedLogMessages.length; i++) {
         partFunc += Math.exp(updatedLogMessages[i]);
         }
         double logPart = Math.log(partFunc);
         for (int i = 0; i < updatedLogMessages.length; i++) {
         updatedLogMessages[i] -= logPart;
         }*/

        return updatedLogMessages;
    }

    double[] subtractReceiver(double[] sumLogMessages, MRFNode sender, MRFNode receiver) {
        double[] sumLogMessagesMinusReceiver = new double[sender.getNumLabels()];
        for (int i = 0; i < sender.getNumLabels(); i++) {
            sumLogMessagesMinusReceiver[i] = sumLogMessages[i] - getLogMessage(receiver, sender, sender.labelList.get(i));
        }
        return sumLogMessagesMinusReceiver;
    }

    double getNonNormalizedLogMessage(double sum, double normalizer) {
        return Math.log(sum) + normalizer;
    }

    void normalizeLogMessages(double[] unNormalized, double largesLogMessage) {
        for (int i = 0; i < unNormalized.length; i++) {
            unNormalized[i] -= largesLogMessage;
        }
    }

    public double sumNormalizedToBeSummed(double[] toBeSummed, double normalizer, int numElements) {
        double sum = 0.;
        for (int sendingLabelIndex = 0; sendingLabelIndex < numElements; sendingLabelIndex++) {
            toBeSummed[sendingLabelIndex] -= normalizer;
            sum += Math.exp(toBeSummed[sendingLabelIndex]);
        }
        return sum;
    }

    public double getMessageUpdateInSum(double pairwisePot, double edgeProbability, double sendingLabelPot, double sumLogMessage) {
        return (((pairwisePot / edgeProbability) + sendingLabelPot) / constRT) + sumLogMessage;
    }

    public int getLabelIndex(MRFNode node, MRFLabel label) {
        return node.labelList.indexOf(label);
    }

    /**
     * Returns a list of the sum log messages from sending to receiving the
     * first entry in the list contains the value of the maximum sum log message
     *
     * @param sendingNode
     * @param receivingNode
     * @return
     */
    double[] getSumLogMessages(MRFNode sendingNode, MRFNode receivingNode) {
        double[] sumLogMessages = new double[sendingNode.labelList.size() + 1];
        //keep maximum in first slot
        int index = 1;
        double max = Double.NEGATIVE_INFINITY;
        for (MRFLabel label : sendingNode.labelList) {
            sumLogMessages[index] = getSumLogMessage(sendingNode, label, receivingNode);
            max = Math.max(sumLogMessages[index], max);
            index++;
        }
        sumLogMessages[0] = max;
        return sumLogMessages;
    }

    /**
     * Computes the log of the righthand side of the equation in (39) in TRBP
     * paper In the non-log-domain this is the product messages raise to the
     * power of the edge probability In log domain this is a sum
     */
    double getSumLogMessage(MRFNode sendingNode, MRFLabel sendingLabel, MRFNode receivingNode
    ) {
        double sum = 0.;
        for (MRFNode nodeV : sendingNode.neighborList) {
            if (!nodeV.equals(receivingNode)) {
                double edgeProbVToSender = getEdgeProbability(sendingNode, nodeV);
                double logMessageVToSender = getLogMessage(nodeV, sendingNode, sendingLabel);
                sum += edgeProbVToSender * logMessageVToSender;
            }
        }
        double edgeProbReceivToSend = getEdgeProbability(sendingNode, receivingNode);
        double logMessageReceivToSend = getLogMessage(receivingNode, sendingNode, sendingLabel);
        sum -= (1 - edgeProbReceivToSend) * logMessageReceivToSend;

        return sum;
    }

    double getSumLogMessageV2(MRFNode receivingNode, MRFLabel label) {
        double sum = 0;
        for (MRFNode nodeV : receivingNode.neighborList) {
            double edgeProbVToSender = getEdgeProbability(receivingNode, nodeV);
            double logMessageVToSender = getLogMessage(nodeV, receivingNode, label);
            sum += edgeProbVToSender * logMessageVToSender;
        }
        return sum;
    }

    double[] getSumLogMessagesV2(MRFNode receivingNode) {
        double[] sumLogMessages = new double[receivingNode.labelList.size()];
        int index = 0;
        for (MRFLabel label : receivingNode.labelList) {
            sumLogMessages[index] = getSumLogMessageV2(receivingNode, label);
            index++;
        }
        return sumLogMessages;
    }

    public int getLabelListSize(MRFNode node) {
        return node.labelList.size();
    }

    public double calcUBLogZ() {
        return -(calcFreeEnergy() + this.emat.getConstTerm()) / this.constRT;
    }

    private double calcFreeEnergy() {
        double enthalpy = getEnthalpy();
        double entropy = getEntropy();
        double freeEnergy = enthalpy - this.constRT * entropy;
        return freeEnergy;
    }

    private double getEnthalpy() {
        double enthalpy = 0.0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode node1 = nodeList.get(i);
            enthalpy += getSingleNodeEnthalpy(node1);
            for (int j = 0; j < i; j++) {
                MRFNode node2 = nodeList.get(j);
                if (this.interactionGraph[node1.index][node2.index]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    public double getEntropy() {
        double entropy = 0.0;
        double[] singleNodeEntropyList = new double[this.nodeList.size()];
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode nodeI = this.nodeList.get(i);
            double singleNodeEntropy = getSingleNodeEntropy(nodeI);
            singleNodeEntropyList[i] = singleNodeEntropy;
            entropy += singleNodeEntropy;
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (interactionGraph[i][j]) {
                    double edgeProb = getEdgeProbability(nodeI, nodeJ);
                    double mutualInf = getMutualInformation(nodeI, nodeJ);
                    double entropyI = singleNodeEntropyList[i];
                    double entropyJ = singleNodeEntropyList[j];
                    //Update Edge Weights to Save Computational Cost
                    this.edgeWeights[i][j] = -mutualInf;///Math.sqrt(entropyI*entropyJ);
                    this.edgeClampWeights[i][j] = -mutualInf;
                    entropy -= edgeProb * mutualInf;
                }
            }
        }
        return entropy;
    }

    private double getSingleNodeEnthalpy(MRFNode node) {
        double enthalpy = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double E = this.emat.getOneBody(node.posNum, node.labelList.get(rot).labelNum);
            double prob = this.marginalProbabilities.getOneBody(node.index, rot);
            enthalpy += E * prob;
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

    private double getPairwiseNodeEnthalpy(MRFNode nodeI, MRFNode nodeJ) {
        double enthalpy = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double E = emat.getPairwise(nodeI.posNum, nodeI.labelList.get(rotI).labelNum, nodeJ.posNum, nodeJ.labelList.get(rotJ).labelNum);
                double prob = this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
                enthalpy += E * prob;
            }
        }
        return enthalpy;
    }

    private void initializeEdgeWeights() {
        this.edgeWeights = new double[this.numNodes][];
        this.edgeClampWeights = new double[this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            this.edgeWeights[i] = new double[i];
            this.edgeClampWeights[i] = new double[i];
        }
    }

    private void updateEdgeProbabilities(int numIter) {
        MinSpanningTree mst = new MinSpanningTree(this.edgeWeights, interactionGraph);
        double[][] descentDirection = mst.mstVector;

        double stepSize = 2.0 / (numIter + 2.0);
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                this.edgeProbabilities[i][j] = stepSize * descentDirection[i][j] + (1 - stepSize) * edgeProbabilities[i][j];
            }
        }
    }

    double getLogMessage(MRFNode sendingNode, MRFNode receivingNode, MRFLabel receivingLabel
    ) {
        return this.logMessages[sendingNode.index][receivingNode.index][receivingNode.labelList.indexOf(receivingLabel)];
    }

    double getPairwisePotential(MRFNode nodeI, MRFLabel labelI, MRFNode nodeJ, MRFLabel labelJ
    ) {
        return -this.emat.getPairwise(nodeI.posNum, labelI.labelNum, nodeJ.posNum, labelJ.labelNum);
    }

    double getPairwisePotential(MRFNode nodeI, MRFLabel labelI, MRFNode nodeJ, MRFLabel labelJ, ArrayList<ArrayList<Double>> energies){
        if (nodeI.posNum > nodeJ.posNum){
            return -energies.get(labelI.labelNum).get(labelJ.labelNum);
        }
        return -energies.get(labelJ.labelNum).get(labelI.labelNum);
    }
    
    double getOneBodyPotential(MRFNode node, MRFLabel label
    ) {
        return -this.emat.getOneBody(node.posNum, label.labelNum);
    }
    
    double getOneBodyPotential(MRFLabel label, ArrayList<Double> energies){
        return -energies.get(label.labelNum);
    }
    
    private double getEdgeProbSum() {
        double sum = 0;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                sum += this.edgeProbabilities[i][j];
            }
        }
        return sum;
    }

    private double getEdgeProbSum(int indexLeaveOut) {
        double sum = 0;
        for (int i = 0; i < this.numNodes; i++) {
            if (i != indexLeaveOut) {
                for (int j = 0; j < i; j++) {
                    if (j != indexLeaveOut) {
                        sum += this.edgeProbabilities[i][j];
                    }
                }
            }
        }
        return sum;
    }

    public double getLogZ() {
        return this.logZ;
    }

    public static void setNumEdgeProbUpdates(int numUpdates) {
        assert (numUpdates >= 0);
        maxNumEdgeUpdates = numUpdates;
    }
}
