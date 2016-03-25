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
public class TreeReweightedBeliefPropagation {

    double logZ;

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
    final int maxIterations = 5000;

    double constRT = PoissonBoltzmannEnergy.constRT;

    double damping = 0.4;

    //p_e are the edge probabilities over spanning trees
    double[][] edgeProbabilities;
    double[][][] messages;

    public ExpFunction ef = new ExpFunction();

    double maxChange;

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

        int numIter = 0;
        for (int i = 0; i < this.maxIterations; i++) {
            numIter++;
            double[][][] messagesNPlus1 = updateMessages(this.messages);
            this.messages = messagesNPlus1;
            if ((maxChange < 1e-8) && (i > 1)) {
                break;
            }
        }
        System.out.println("TRBP took: " + numIter + " iterations");
        updateMarginals(this.messages);
        checkMarginals();
    }

    private double[][][] updateMessages(double[][][] previousMessages) {
        double[][][] updatedMessages = CreateMatrix.create3DMsgMat(numNodes, numLabelsPerNode, 0.0);
        double maxChangeMessage = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {

                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                if (this.interactionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        double messageAtRotI = 0.0;
                        double expNormI = Double.NEGATIVE_INFINITY;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
                            double toBeExponentiate = (pairwiseE / edgeProbability) + rotJE;
                            expNormI = Math.max(expNormI, toBeExponentiate);
                        }
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.nodeNum, rotJ.labelNum);
                            double normalized = ((pairwiseE / edgeProbability) + rotJE) - expNormI;
//                            double messageFromRotJ = Math.exp(normalized/this.constRT);
                            double messageFromRotJ = Math.exp(normalized);
                            messageFromRotJ = messageFromRotJ * (getProductMessages(nodeJ, rotJindex, nodeI, previousMessages));
                            messageAtRotI = messageAtRotI + (messageFromRotJ);
                        }
                        updatedMessages[j][i][rotIindex] = messageAtRotI;
                        partFunctionI += messageAtRotI;
                    }
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        updatedMessages[j][i][rotIindex] /= partFunctionI;
                        updatedMessages[j][i][rotIindex] = damping * updatedMessages[j][i][rotIindex] + ((1 - damping) * previousMessages[j][i][rotIindex]);
                        double prevMessage = previousMessages[j][i][rotIindex];
                        double change = Math.abs(prevMessage - updatedMessages[j][i][rotIindex]);
                        maxChangeMessage = Math.max(change, maxChangeMessage);
                    }

                    //Now we send message from i->j
                    double partFunctionJ = 0.0;
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        double messageAtRotJ = 0.0;
                        double expNormJ = Double.NEGATIVE_INFINITY;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);

                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotIE = -this.emat.getOneBody(nodeI.nodeNum, rotI.labelNum);

                            double toBeExponentiated = (pairwiseE / edgeProbability) + rotIE;
                            expNormJ = Math.max(expNormJ, toBeExponentiated);
                        }
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);

                            double pairwiseE = -this.emat.getPairwise(nodeI.nodeNum, rotI.labelNum, nodeJ.nodeNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotIE = -this.emat.getOneBody(nodeI.nodeNum, rotI.labelNum);
                            double normalized = ((pairwiseE / edgeProbability) + rotIE) - expNormJ;
//                            double messageFromRotI = Math.exp(normalized/this.constRT);
                            double messageFromRotI = Math.exp(normalized);
                            messageFromRotI = messageFromRotI * (getProductMessages(nodeI, rotIindex, nodeJ, previousMessages));
                            messageAtRotJ += messageFromRotI;
                        }
                        updatedMessages[i][j][rotJindex] = messageAtRotJ;
                        partFunctionJ += messageAtRotJ;
                    }
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        updatedMessages[i][j][rotJindex] /= partFunctionJ;
                        updatedMessages[i][j][rotJindex] = damping * updatedMessages[i][j][rotJindex] + ((1 - damping) * previousMessages[i][j][rotJindex]);
                        double prevMessage = previousMessages[i][j][rotJindex];
                        double change = Math.abs(prevMessage - updatedMessages[i][j][rotJindex]);
                        maxChangeMessage = Math.max(change, maxChangeMessage);
                    }
                }
            }
        }
        this.maxChange = maxChangeMessage;
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
        double partFuncOverall = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            double partFunc = 0.0;
            double expNorm = Double.NEGATIVE_INFINITY;
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double pairwiseE = this.emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double edgeProb = getEdgeProbability(nodeI.nodeNum, nodeJ.nodeNum);
                double nodeIE = this.emat.getOneBody(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum);
                double nodeJE = this.emat.getOneBody(nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);

                double toBeExponentiated = -((pairwiseE / edgeProb) + nodeIE + nodeJE);
                expNorm = Math.max(expNorm, toBeExponentiated);
            }
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double pairwiseE = this.emat.getPairwise(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum, nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);
                double edgeProb = getEdgeProbability(nodeI.nodeNum, nodeJ.nodeNum);
                double nodeIE = this.emat.getOneBody(nodeI.nodeNum, nodeI.labelList.get(rotI).labelNum);
                double nodeJE = this.emat.getOneBody(nodeJ.nodeNum, nodeJ.labelList.get(rotJ).labelNum);

                double normalized = (-((pairwiseE / edgeProb) + nodeIE + nodeJE)) - expNorm;
//                double marginal = Math.exp(normalized/this.constRT);
                double marginal = Math.exp(normalized);
                marginal = marginal * (getProductMessages(nodeI, rotI, nodeJ, messages));
                marginal = marginal * (getProductMessages(nodeJ, rotJ, nodeI, messages));

                partFunc += marginal;
                this.marginalProbabilies.setPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ, marginal);
            }

            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double unNormalized = this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
                double marginalAtRotI = this.marginalProbabilies.getOneBody(nodeI.nodeNum, rotI);
                double normalized;
                normalized = unNormalized * marginalAtRotI / partFunc;
                partFuncOverall += normalized;
                this.marginalProbabilies.setPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ, normalized);
            }

        }
/*        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double unNormalized = this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
                double normalized = unNormalized / partFuncOverall;
                this.marginalProbabilies.setPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ, normalized);
            }
        } */
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
            sum += this.marginalProbabilies.getOneBody(node.nodeNum, rot);
        }
        if (Math.abs(sum - 1) > 1e-6) {
            throw new RuntimeException("Marginal Not Normalized at Node: " + node.nodeNum);
        }
    }

    private void checkPairwiseMarginal(MRFNode nodeI, MRFNode nodeJ) {
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            double marginal = this.marginalProbabilies.getOneBody(nodeI.nodeNum, rotI);
            double sum = 0.0;
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                sum += this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
            }
            if (Math.abs(sum - marginal) > 1e-6) {
                throw new RuntimeException("Marginal Not Normalized at Node Pairs: " + nodeI.nodeNum + " " + nodeJ.nodeNum);
            }
        }
        for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
            double marginal = this.marginalProbabilies.getOneBody(nodeJ.nodeNum, rotJ);
            double sum = 0.0;
            for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                sum += this.marginalProbabilies.getPairwise(nodeI.nodeNum, rotI, nodeJ.nodeNum, rotJ);
            }
            if (Math.abs(sum - marginal) > 1e-6) {
                throw new RuntimeException("Marginal Not Normalized at Node Pairs: " + nodeJ.nodeNum + " " + nodeI.nodeNum);
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
//            double update = Math.exp(normalized/this.constRT);
            double update = Math.exp(normalized);
            for (MRFNode neighbor : node.neighborList) {
                double message = this.messages[neighbor.nodeNum][node.nodeNum][rot];
                double prob = getEdgeProbability(neighbor.nodeNum, node.nodeNum);
                update = update * (Math.pow(message, prob));
            }
            this.marginalProbabilies.setOneBody(node.nodeNum, rot, update);
            partFunc += update;
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double unNormalized = this.marginalProbabilies.getOneBody(node.nodeNum, rot);
            this.marginalProbabilies.setOneBody(node.nodeNum, rot, unNormalized / partFunc);
        }
    }

    private double getProductMessages(MRFNode nodeReceiving, int labelIndex, MRFNode nodeExcluded, double[][][] currentMessages) {
        double product = 1;
        for (MRFNode node : nodeReceiving.neighborList) {
            if (!(node.equals(nodeExcluded))) {
                double message = currentMessages[node.nodeNum][nodeReceiving.nodeNum][labelIndex];
                double probability = getEdgeProbability(nodeReceiving.nodeNum, node.nodeNum);
                product = product * (Math.pow(message, probability));
            }
        }
        double messageNodeExcluded = currentMessages[nodeExcluded.nodeNum][nodeReceiving.nodeNum][labelIndex];
        double probability = getEdgeProbability(nodeExcluded.nodeNum, nodeReceiving.nodeNum);
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
//                    edgesFromI[j] = 1.0 / ((double) numEdges);
//                    edgesFromI[j] = 0.5;
                    edgesFromI[j] = 2.0/((double) this.nodeList.size());
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
            entropy += (-1.0) * prob * Math.log(prob);
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

    public double getEntropy() {
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

    public double calcUBLogZ() {
        return -(calcFreeEnergy() + this.emat.getConstTerm()) / this.constRT;
    }
}
