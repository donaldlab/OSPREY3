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

    double[][] expNormMessages;
    double[][] expNormMarginals;
    boolean dampLog = false;

    public ExpFunction ef = new ExpFunction();

    double maxChange;
    double averageChange;

    boolean verbose;

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

        this.edgeProbabilities = initializeEdgeProbabilities(this.nonClampledInteractionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);

        runTRBPSeq();
    }

    double[][][] initializeLogMessages(double initVal) {
        return CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, initVal);
    }

    double[][] initializeEdgeProbabilities(boolean[][] interactionGraph) {
        return GraphUtils.getEdgeProbabilities(interactionGraph);
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

    private double getEdgeProbability(MRFNode nodeI, MRFNode nodeJ) {
        int indexNodeI = nodeI.index;
        int indexNodeJ = nodeJ.index;
        if (indexNodeI < indexNodeJ) {
            return this.edgeProbabilities[indexNodeJ][indexNodeI];
        }
        return this.edgeProbabilities[indexNodeI][indexNodeJ];
    }

    private double[][][] updateMessagesSequentially(double[][][] previousLogMessages) {
        double[][][] updatedLogMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);
        double averageChangeInMessage = 0;
        int numMessagesUpdated = 0;

        //Run in sequential order forward then backword
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                if (this.interactionGraph[i][j]) {
                    MRFNode nodeI = this.nodeList.get(i);
                    MRFNode nodeJ = this.nodeList.get(j);
                    double[] logMessagesNodeJToNodeI = getUpdatedLogMessage(nodeJ, nodeI);
                }
            }
        }

    }

    double[] getUpdatedLogMessage(MRFNode sendingNode, MRFNode receivingNode) {
        double[] updatedLogMessages = new double[receivingNode.labelList.size()];
        double largestLogMessage = Double.NEGATIVE_INFINITY;
        for (MRFLabel receivingLabel : receivingNode.labelList) {
            double sum = 0.;
            for (MRFLabel sendingLabel : sendingNode.labelList) {
                //compute everything between the curly brackets in eq 39. of TRBP Paper
                //We do this in the log domain and exponentiate after
                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel);
                double edgeProbability = getEdgeProbability(receivingNode, sendingNode);
                double sendingLabelPot = getOneBodyPotential(sendingNode, sendingLabel);

                double normalized = ((pairwisePot / edgeProbability) + sendingLabelPot - getExpNormMessages(sendingNode, receivingNode));
                double sumLogMessages = getSumLogMessage(sendingNode, sendingLabel, receivingNode);

                double exponential = Math.exp(normalized + sumLogMessages);

                sum += exponential;
            }
            double nonNormalizedLogMessage = Math.log(sum);
            largestLogMessage = Math.max(largestLogMessage, nonNormalizedLogMessage);
            int messageIndex = receivingNode.labelList.indexOf(receivingLabel);
            updatedLogMessages[messageIndex] = nonNormalizedLogMessage;
        }
        //Normalize by subtracting the largest logMessage;
        for (int i = 0; i < updatedLogMessages.length; i++) {
            updatedLogMessages[i] -= largestLogMessage;
        }
        return updatedLogMessages;
    }

    /**
     * Computes the log of the righthand side of the equation in (39) in TRBP
     * paper In the non-log-domain this is the product messages raise to the
     * power of the edge probability In log domain this is a sum
     */
    double getSumLogMessage(MRFNode sendingNode, MRFLabel sendingLabel, MRFNode receivingNode) {
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

    double getLogMessage(MRFNode sendingNode, MRFNode receivingNode, MRFLabel receivingLabel) {
        return this.logMessages[sendingNode.index][receivingNode.index][receivingNode.labelList.indexOf(receivingLabel)];
    }

    double getExpNormMessages(MRFNode sendingNode, MRFNode receivingNode) {
        return this.expNormMessages[sendingNode.index][receivingNode.index];
    }

    double getPairwisePotential(MRFNode nodeI, MRFLabel labelI, MRFNode nodeJ, MRFLabel labelJ) {
        return -this.emat.getPairwise(nodeI.nodeNum, labelI.labelNum, nodeJ.nodeNum, labelJ.labelNum);
    }

    double getOneBodyPotential(MRFNode node, MRFLabel label) {
        return -this.emat.getOneBody(node.nodeNum, label.labelNum);
    }
}
