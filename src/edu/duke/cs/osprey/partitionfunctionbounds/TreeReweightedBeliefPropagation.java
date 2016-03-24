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
import java.math.BigDecimal;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class TreeReweightedBeliefPropagation {

    BigDecimal partitionFunction;

    //Node list
    ArrayList<MRFNode> nodeList;
    //Emat defines the potential functions
    int numNodes;
    int[] numLabelsPerNode;

    EnergyMatrix emat;
    //interactionGraph;
    boolean[][] interactionGraph;

    TupleMatrix<Double> marginalProbabilies;

    //threshold for convergence
    double threshold = 1e-8;
    final int maxIterations = 1000;

    double constRT = PoissonBoltzmannEnergy.constRT;

    //p_e are the edge probabilities over spanning trees
    double[][] edgeProbabilities;
    double[][][] messages;

    public ExpFunction ef = new ExpFunction();

    public TreeReweightedBeliefPropagation(MarkovRandomField mrf) {
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
        this.messages = initializeMessages(1.0);

        this.marginalProbabilies = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);

        for (int i = 0; i < this.maxIterations; i++) {
            double[][][] messagesNPlus1 = updateMessages(this.messages);
            this.messages = messagesNPlus1;
            if (i > 50) {
                updateMarginals(this.messages);
                double freeEnergy = calcFreeEnergy();
            }
        }
    }

    private double[][][] updateMessages(double[][][] previousMessages) {
        double[][][] updatedMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {

                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                if (this.interactionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        BigDecimal messageAtRotI = new BigDecimal(0.0);
                        double expNorm = Double.NEGATIVE_INFINITY;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
                            double toBeExponentiate = (pairwiseE/edgeProbability) + rotJE;
                            expNorm = Math.max(expNorm,toBeExponentiate);
                        }
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
                            double normalize = ((pairwiseE/edgeProbability)+rotJE) - expNorm;
                            double messageFromRotJ = Math.pow(Math.E,normalize);
                            messageFromRotJ = messageFromRotJ.multiply(getProductMessages(nodeJ, rotJindex, nodeI, previousMessages));
                            messageAtRotI = messageAtRotI.add(messageFromRotJ);
                        }
                        double messageAtRotId = messageAtRotI.doubleValue();
                        updatedMessages[j][i][rotIindex] = messageAtRotId;
                        partFunctionI += messageAtRotId;
                    }
                    double partFunc2 = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        updatedMessages[j][i][rotIindex] /= partFunctionI;
                        if (updatedMessages[j][i][rotIindex] < 1e-6) {
                            updatedMessages[j][i][rotIindex] = 1e-6;
                        }
                        partFunc2 += updatedMessages[j][i][rotIindex];
                    }
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        updatedMessages[j][i][rotIindex] /= partFunc2;
                    }

                    //Now we send message from i->j
                    double partFunctionJ = 0.0;
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        BigDecimal messageAtRotJ = new BigDecimal(0.0);
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);

                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotIE = -this.emat.getOneBody(nodeI.nodeNum, rotI.labelNum);
                            BigDecimal messageFromRotI = this.ef.exp((pairwiseE / edgeProbability) + rotIE);
                            messageFromRotI = messageFromRotI.multiply(getProductMessages(nodeI, rotIindex, nodeJ, previousMessages));
                            messageAtRotJ = messageAtRotJ.add(messageFromRotI);
                        }
                        double messageAtRotJd = messageAtRotJ.doubleValue();
                        updatedMessages[i][j][rotJindex] = messageAtRotJd;
                        partFunctionJ += messageAtRotJd;
                    }
                    double partFuncJ2 = 0.0;
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        updatedMessages[i][j][rotJindex] /= partFunctionJ;
                        if (updatedMessages[i][j][rotJindex] < 1e-6) {
                            updatedMessages[i][j][rotJindex] = 1e-6;
                        }
                        partFuncJ2 += updatedMessages[i][j][rotJindex];
                    }
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        updatedMessages[i][j][rotJindex] /= partFuncJ2;
                    }
                }
            }
        }
        return updatedMessages;
    }

    private void updateMarginals(double[][][] messages) {
        //first update node marginals
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            updateNodeMarginal(node, messages);
        }
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                updateEdgeMarginal(nodeI, nodeJ, messages);
            }
        }
    }

    private void updateEdgeMarginal(MRFNode nodeI, MRFNode nodeJ, double[][][] messages) {
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            double partFunc = 0.0;
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double pairwiseE = this.emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double edgeProb = getEdgeProbability(nodeI.nodeNum, nodeJ.nodeNum);
                double nodeIE = this.emat.getOneBody(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum);
                double nodeJE = this.emat.getOneBody(nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);

                BigDecimal marginal = this.ef.exp(-((pairwiseE / edgeProb) + nodeIE + nodeJE));
                marginal = marginal.multiply(getProductMessages(nodeI, rotI, nodeJ, messages));
                marginal = marginal.multiply(getProductMessages(nodeJ, rotJ, nodeI, messages));

                double marginald = marginal.doubleValue();
                double marginalAtRotI = this.marginalProbabilies.getOneBody(nodeI.nodeNum, rotI);
                double marginalAtRotJ = this.marginalProbabilies.getOneBody(nodeJ.nodeNum, rotJ);
                if ((marginalAtRotI == 0) || (marginalAtRotJ == 0)) {
                    marginald = 0.0;
                }
                partFunc += marginald;
                this.marginalProbabilies.setPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ, marginald);
            }
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double unNormalized = this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
                double marginalAtRotI = this.marginalProbabilies.getOneBody(nodeI.nodeNum, rotI);
                double marginalAtRotJ = this.marginalProbabilies.getOneBody(nodeJ.nodeNum, rotJ);
                double normalized;
                if (partFunc == 0.0) {
                    if ((marginalAtRotI == 0) || (marginalAtRotJ == 0)) {
                        normalized = 0.0;
                    } else {
                        normalized = marginalAtRotI / nodeJ.labelList.size();
                    }
                } else {
                    normalized = unNormalized * marginalAtRotI / partFunc;
                }
                this.marginalProbabilies.setPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ, normalized);
            }
        }
    }

    private void updateNodeMarginal(MRFNode node, double[][][] messages) {
        double partFunc = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            BigDecimal update = this.ef.exp(-this.emat.getOneBody(node.nodeNum, node.labelList.get(rot).labelNum));
            for (MRFNode neighbor : node.neighborList) {
                double message = this.messages[neighbor.nodeNum][node.nodeNum][rot];
                double prob = getEdgeProbability(neighbor.nodeNum, node.nodeNum);
                update = update.multiply(new BigDecimal(Math.pow(message, prob)));
            }
            double updateD = update.doubleValue();
            this.marginalProbabilies.setOneBody(node.nodeNum, rot, updateD);
            partFunc += updateD;
        }
        double partFunc2 = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double unNormalized = this.marginalProbabilies.getOneBody(node.nodeNum, rot);
            double normalized = unNormalized / partFunc;
            if (normalized < 1e-6) {
                normalized = 1e-6;
            }
            partFunc2 += normalized;
            this.marginalProbabilies.setOneBody(node.nodeNum, rot, unNormalized / partFunc);
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double unNormalized = this.marginalProbabilies.getOneBody(node.nodeNum, rot);
            this.marginalProbabilies.setOneBody(node.nodeNum, rot, unNormalized / partFunc2);
        }
    }

    private BigDecimal getProductMessages(MRFNode nodeReceiving, int labelIndex, MRFNode nodeExcluded, double[][][] currentMessages) {
        BigDecimal product = new BigDecimal("1.0");
        for (MRFNode node : nodeReceiving.neighborList) {
            if (!(node.equals(nodeExcluded))) {
                double message = currentMessages[node.nodeNum][nodeReceiving.nodeNum][labelIndex];
                double probability = getEdgeProbability(nodeReceiving.nodeNum, node.nodeNum);
                product = product.multiply(new BigDecimal(Math.pow(message, probability)));
            }
        }
        double messageNodeExcluded = currentMessages[nodeExcluded.nodeNum][nodeReceiving.nodeNum][labelIndex];
        double probability = getEdgeProbability(nodeExcluded.nodeNum, nodeReceiving.nodeNum);
        double bottom = Math.pow(messageNodeExcluded, (1 - probability));
        if (bottom == 0.0) {
            product = product;
        } else {
            product = product.divide(new BigDecimal(Math.pow(messageNodeExcluded, (1 - probability))), this.ef.mc);
        }
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

        return CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, initVal);
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
                    edgesFromI[j] = 1.0 / ((double) numEdges);
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
            double prob = this.marginalProbabilies.getOneBody(node.nodeNum, rot);
            enthalpy += E * prob;
        }
        return enthalpy;
    }

    private double getPairwiseNodeEnthalpy(MRFNode nodeI, MRFNode nodeJ) {
        double enthalpy = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double E = emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double prob = this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
                if (Double.isNaN(E * prob)) {
                    enthalpy += E * prob;
                }
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
                if (this.interactionGraph[node1.nodeNum][node2.nodeNum]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double prob = this.marginalProbabilies.getOneBody(node.nodeNum, rot);
            if (prob == 0.0) {
                entropy += 0.0;
            } else {
                entropy += (-1.0) * prob * Math.log(prob);
            }
        }
        return entropy;
    }

    private double getMutualInformation(MRFNode nodeI, MRFNode nodeJ) {
        double mutualInf = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double probIJ = this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
                double probI = this.marginalProbabilies.getOneBody(nodeI.nodeNum, rotI);
                double probJ = this.marginalProbabilies.getOneBody(nodeJ.nodeNum, rotJ);
                if (!(probIJ == 0.0)) {
                    mutualInf += probIJ * Math.log(probIJ / (probI * probJ));
                }
            }
        }
        return mutualInf;
    }

    private double getEntropy() {
        double entropy = 0.0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            MRFNode nodeI = this.nodeList.get(i);
            entropy += getSingleNodeEntropy(nodeI);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (interactionGraph[i][j]) {
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
}
