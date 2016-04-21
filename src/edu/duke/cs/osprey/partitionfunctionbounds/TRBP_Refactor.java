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
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class TRBP_Refactor {

    double logZ = Double.POSITIVE_INFINITY;

    ArrayList<MRFNode> nodeList;
    int numNodes;
    int[] numLabelsPerNode;

    UpdatedEmat emat;
    boolean[][] interactionGraph;

    TupleMatrix<Double> marginalProbabilities;

    double threshold = 1e-6;
    double constRT = PoissonBoltzmannEnergy.constRT;

    double damping = 0.5;
    int numEdgeProbUpdates;

    public double[][] edgeProbabilities;
    double[][][] logMessages;
//    double[][][] messages;

    double[][][] expNormMessagesAtRot;
    double[][] expNormMessages;
    double[][] expNormPairMarginals;
    double[] expNormNodeMarginals;

    boolean useLogDomain = true;
    boolean debug = true;

    public ExpFunction ef = new ExpFunction();

    double maxChange;
    double averageChange;

    boolean verbose = true;

    public TRBP_Refactor(ReparamMRF mrf) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.nonClampedInteractionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.labelList.size();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.interactionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);

        runTRBP();
    }

    private void runTRBP() {

        computeExpNormals();
        computeExpNormals2();
        for (int i = 0; i < 40; i++) {
            updateMessagesSequentially(logMessages);
            updateMarginals();
            double currentlogZ = calcUBLogZ();
            if (verbose) {
                System.out.println("Average Change in Messages: " + this.maxChange + "  LogZUB: " + currentlogZ);
            }
        }
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

    private void computeExpNormals2() {
        this.expNormMessagesAtRot = new double[this.numNodes][this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (this.interactionGraph[i][j]) {
                    MRFNode nodeI = this.nodeList.get(i);
                    MRFNode nodeJ = this.nodeList.get(j);
                    expNormMessagesAtRot[i][j] = getExpNormMessage(nodeI, nodeJ);
                    expNormMessagesAtRot[j][i] = getExpNormMessage(nodeJ, nodeI);
                }
            }
        }

    }

    private double[] getExpNormMessage(MRFNode sendingNode, MRFNode receivingNode) {
        double[] expNorm = new double[receivingNode.labelList.size()];
        for (int rot = 0; rot < expNorm.length; rot++) {
            MRFLabel label = receivingNode.labelList.get(rot);
            double maxPotential = Double.NEGATIVE_INFINITY;
            for (int rotS = 0; rotS < sendingNode.labelList.size(); rotS++) {
                MRFLabel labelS = sendingNode.labelList.get(rotS);
                double pairPot = getPairwisePotential(sendingNode, labelS, receivingNode, label);
                double edgeProb = getEdgeProbability(sendingNode, receivingNode);
                double labelSPot = getOneBodyPotential(sendingNode, labelS);
                maxPotential = Math.max(maxPotential, (pairPot / edgeProb) + labelSPot);
            }
            expNorm[rot] = maxPotential;
        }
        return expNorm;
    }

    private void computeExpNormals() {
        double[][] expNormMessage = new double[this.numNodes][this.numNodes];
        double[][] expNormPairMarginal = new double[this.numNodes][this.numNodes];
        double[] expNormNodeMarginal = new double[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode nodeI = this.nodeList.get(i);
            for (int nodeMarginal = 0; nodeMarginal < nodeI.labelList.size(); nodeMarginal++) {
                MRFLabel label = nodeI.labelList.get(nodeMarginal);
                double labelPot = getOneBodyPotential(nodeI, label);
                expNormNodeMarginal[i] = Math.max(expNormNodeMarginal[i], labelPot);
            }
        }

        for (int i = 0; i < this.numNodes; i++) {
            MRFNode nodeI = this.nodeList.get(i);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        MRFLabel rotI = nodeI.labelList.get(rotIindex);
                        MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                        double pairPot = getPairwisePotential(nodeI, rotI, nodeJ, rotJ);
                        double edgeProbability = getEdgeProbability(nodeI, nodeJ);
                        double labelJPot = getOneBodyPotential(nodeJ, rotJ);
                        double labelIPot = getOneBodyPotential(nodeI, rotI);
                        double toBeExponentiateJI = (pairPot / edgeProbability) + labelJPot;
                        double toBeExponentiateIJ = (pairPot / edgeProbability) + labelIPot;
                        double toBeExponentiatedMarginal = (pairPot / edgeProbability) + labelIPot + labelJPot;
                        expNormMessage[j][i] = Math.max(expNormMessage[j][i], toBeExponentiateJI);
                        expNormMessage[i][j] = Math.max(expNormMessage[i][j], toBeExponentiateIJ);
                        expNormPairMarginal[i][j] = Math.max(expNormPairMarginal[i][j], toBeExponentiatedMarginal);
                        expNormPairMarginal[j][i] = expNormPairMarginal[i][j];
                    }
                }
            }
        }
        this.expNormMessages = expNormMessage;
        this.expNormPairMarginals = expNormPairMarginal;
        this.expNormNodeMarginals = expNormNodeMarginal;
    }

    private double getEdgeProbability(MRFNode nodeI, MRFNode nodeJ) {
        int indexNodeI = nodeI.index;
        int indexNodeJ = nodeJ.index;
        if (indexNodeI < indexNodeJ) {
            return this.edgeProbabilities[indexNodeJ][indexNodeI];
        }
        return this.edgeProbabilities[indexNodeI][indexNodeJ];
    }

    private void updateMessagesSequentially(double[][][] previousMessages) {
        double averageChangeInMessage = 0;
        int numMessagesUpdated = 0;

        int[][] messagesOrdering = getMessagePassingOrdering();

        for (int[] nodePair : messagesOrdering) {
            //Node pair consists of two ints. The first indexes to the sending node
            //The second indexes to the receiving node
            MRFNode nodeI = this.nodeList.get(nodePair[0]);
            MRFNode nodeJ = this.nodeList.get(nodePair[1]);

            //Get updated log Messages from nodeI to nodeJ
            double[] logMessagesNodeIToNodeJ = getUpdatedLogMessage(nodeI, nodeJ);
            double[] previousLogMessages = previousMessages[nodeI.index][nodeJ.index];
            //double[] messagesNodeJToNodeI = getUpdatedMessages(nodeJ, nodeI);

            //Damp Messages in Log Domain and keep track of change in messages
            for (int logMessage = 0; logMessage < logMessagesNodeIToNodeJ.length; logMessage++) {
                //Get previous log message to keep track of change
                double previousLogMessage = previousLogMessages[logMessage];
                //Dampen messages to get the new log message
                previousLogMessages[logMessage] = damping * logMessagesNodeIToNodeJ[logMessage]
                        + (1 - damping) * previousLogMessage;
                //Get new log message
                double newLogMessage = previousLogMessages[logMessage];
                //Get change in log messages
                double change = Math.abs(newLogMessage - previousLogMessage);
                averageChangeInMessage += change;
                numMessagesUpdated++;
            }
        }
        //Update average change
        this.averageChange = averageChangeInMessage / numMessagesUpdated;
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

            double labelPot = getOneBodyPotential(node, label)/this.constRT;
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
            partFunc += Math.exp(marginals[labelIndex]);
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
        for (int labelIndexI = 0; labelIndexI < nodeI.labelList.size(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.labelList.size(); labelIndexJ++) {
                MRFLabel labelI = nodeI.labelList.get(labelIndexI);
                MRFLabel labelJ = nodeJ.labelList.get(labelIndexJ);
                double pairPot = getPairwisePotential(nodeI, labelI, nodeJ, labelJ);
                double edgeProb = getEdgeProbability(nodeI, nodeJ);
                double labelIPot = getOneBodyPotential(nodeI, labelI);
                double labelJPot = getOneBodyPotential(nodeJ, labelJ);
                
                double potential = ((pairPot / edgeProb) + labelIPot + labelJPot)/this.constRT;
//                double normalizedPot = ((pairPot / edgeProb) + labelIPot + labelJPot - this.expNormPairMarginals[nodeI.index][nodeJ.index]) / this.constRT;
                double sumWeightedLogMessages = getSumLogMessage(nodeI, labelI, nodeJ) + getSumLogMessage(nodeJ, labelJ, nodeI);
                normalizer = Math.max(normalizer, potential + sumWeightedLogMessages);
                
//                double nonNormalizedLogMarginal = normalizedPot + sumWeightedLogMessages;
//                double nonNormalizedMarginal = Math.exp(nonNormalizedLogMarginal);
//                partFunc += nonNormalizedMarginal;
                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, potential + sumWeightedLogMessages);
            }
        }
        for (int labelIndexI = 0; labelIndexI < nodeI.labelList.size(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.labelList.size(); labelIndexJ++) {
                //SUBTRACT OUT NORMALIZER!!!!!!!!!!!!!!!!
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
    int[][] getMessagePassingOrdering() {
        int numMessages = 2 * getNumEdges(this.interactionGraph);
        int[][] messagePassingOrdering = new int[numMessages][];

        int currentMessage = 0;
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
        return messagePassingOrdering;
    }

    double[] getUpdatedLogMessage(MRFNode sendingNode, MRFNode receivingNode
    ) {
        double[] updatedLogMessages = new double[receivingNode.labelList.size()];
        double largestLogMessage = Double.NEGATIVE_INFINITY;
        for (MRFLabel receivingLabel : receivingNode.labelList) {
            int index = receivingNode.labelList.indexOf(receivingLabel);
            double normalizer = Double.NEGATIVE_INFINITY;
            double[] toBeSummed = new double[sendingNode.labelList.size()];
            for (int sendingLabelIndex = 0; sendingLabelIndex < sendingNode.labelList.size(); sendingLabelIndex++) {
                MRFLabel sendingLabel = sendingNode.labelList.get(sendingLabelIndex);

                //compute everything between the curly brackets in eq 39. of TRBP Paper
                //We do this in the log domain and exponentiate aftern
                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel);
                double edgeProbability = getEdgeProbability(receivingNode, sendingNode);
                double sendingLabelPot = getOneBodyPotential(sendingNode, sendingLabel);

//                double normalized = ((pairwisePot / edgeProbability) + sendingLabelPot - getExpNormMessagesAtRot(sendingNode, receivingNode, index));
                double sumLogMessages = getSumLogMessage(sendingNode, sendingLabel, receivingNode);
                toBeSummed[sendingLabelIndex] = (((pairwisePot / edgeProbability) + sendingLabelPot) / constRT) + sumLogMessages;
                normalizer = Math.max(normalizer, toBeSummed[sendingLabelIndex]);
                //              double exponential = Math.exp((normalized / constRT) + sumLogMessages);
                //                sum += exponential;
            }
            double sum = 0.;
            for (int sendingLabelIndex = 0; sendingLabelIndex < sendingNode.labelList.size(); sendingLabelIndex++) {
                toBeSummed[sendingLabelIndex] -= normalizer;
                sum += Math.exp(toBeSummed[sendingLabelIndex]);
            }
            double nonNormalizedLogMessage = Math.log(sum) + normalizer;
            largestLogMessage = Math.max(largestLogMessage, nonNormalizedLogMessage);
            int messageIndex = receivingNode.labelList.indexOf(receivingLabel);
            updatedLogMessages[messageIndex] = nonNormalizedLogMessage;
        }
        //Normalize by subtracting the largest logMessage;
        for (int i = 0; i < updatedLogMessages.length; i++) {
            updatedLogMessages[i] -= largestLogMessage;
        }
        double partFunc = 0.0;
        for (int i = 0; i < updatedLogMessages.length; i++) {
            partFunc += Math.exp(updatedLogMessages[i]);
        }
        double logPart = Math.log(partFunc);
        for (int i = 0; i < updatedLogMessages.length; i++) {
            updatedLogMessages[i] -= logPart;
        }
        double sum = 0;
        for (double logMessage : updatedLogMessages) {
            sum += Math.exp(logMessage);
        }
        return updatedLogMessages;
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
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode nodeI = this.nodeList.get(i);
            entropy += getSingleNodeEntropy(nodeI);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (interactionGraph[i][j]) {
                    double edgeProb = getEdgeProbability(nodeI, nodeJ);
                    entropy -= edgeProb * getMutualInformation(nodeI, nodeJ);
                }
            }
        }
        return entropy;
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
                double E = emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double prob = this.marginalProbabilities.getPairwise(nodeI.index, rotI, nodeJ.index, rotJ);
                enthalpy += E * prob;
            }
        }
        return enthalpy;
    }

    double getLogMessage(MRFNode sendingNode, MRFNode receivingNode, MRFLabel receivingLabel
    ) {
        return this.logMessages[sendingNode.index][receivingNode.index][receivingNode.labelList.indexOf(receivingLabel)];
    }

    double getExpNormMessages(MRFNode sendingNode, MRFNode receivingNode
    ) {
        return this.expNormMessages[sendingNode.index][receivingNode.index];
    }

    double getExpNormMessagesAtRot(MRFNode sendingNode, MRFNode receivingNode, int receivingLabelIndex
    ) {
        return this.expNormMessagesAtRot[sendingNode.index][receivingNode.index][receivingLabelIndex];
    }

    double getPairwisePotential(MRFNode nodeI, MRFLabel labelI, MRFNode nodeJ, MRFLabel labelJ
    ) {
        return -this.emat.getPairwise(nodeI.nodeNum, labelI.labelNum, nodeJ.nodeNum, labelJ.labelNum);
    }

    double getOneBodyPotential(MRFNode node, MRFLabel label
    ) {
        return -this.emat.getOneBody(node.nodeNum, label.labelNum);
    }
}
