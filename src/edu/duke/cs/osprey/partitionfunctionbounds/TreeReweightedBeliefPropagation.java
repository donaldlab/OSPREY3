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

    double logZ = Double.POSITIVE_INFINITY;

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
    final int maxIterations = 10000;

    double constRT = PoissonBoltzmannEnergy.constRT;

    double damping = 0.4;

    //p_e are the edge probabilities over spanning trees
    public double[][] edgeProbabilities;
    double[][][] messages;

    int numEdgeProbUpdates = 20;

    double[][] expNormMessages;
    double[][] expNormMarginals;
    boolean dampLog = false;

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

        runTRBP();
    }

    private void runTRBP() {
        for (int j = 0; j < numEdgeProbUpdates; j++) {
            computeExpNormals();
            int numIter = 0;
            while (maxChange > threshold || numIter < 2) {
                maxChange = 0;
                numIter++;
                double[][][] messagesNPlus1 = updateMessages(this.messages);
                checkMessages(messagesNPlus1);
                this.messages = messagesNPlus1;
            }
//            System.out.println("TRBP took: " + numIter + " iterations");
            updateMarginals(this.messages);
            checkMarginals();
//            System.out.println("LogZUB: "+logZ);

            double currentlogZ = calcUBLogZ();
//            System.out.println("LogZUB: " + currentlogZ);
            if (!Double.isNaN(currentlogZ)) {
                this.logZ = Math.min(this.logZ, currentlogZ);
            }
            if (j < numEdgeProbUpdates - 1) {
                MinSpanningTree mst = new MinSpanningTree(getEdgeWeights(), interactionGraph);
                updateEdgeProbabilies(mst.mstVector, j + 1);
            }
        }
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
                        double pairwiseE = -this.emat.getPairwise(nodeI.posNum, rotI.labelNum, nodeJ.posNum, rotJ.labelNum);
                        double edgeProbability = getEdgeProbability(i, j);
                        double rotJE = -this.emat.getOneBody(nodeJ.posNum, rotJ.labelNum);
                        double rotIE = -this.emat.getOneBody(nodeI.posNum, rotI.labelNum);
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
                if (this.interactionGraph[i][j]) {
                    //First we send message from j->i
                    double partFunctionI = 0.0;
                    for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                        double messageAtRotI = 0.0;
                        for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);
                            //nodeI.nodeNum should be the same as i, nodeJ.nodeNum should be the same as j
                            double pairwiseE = -this.emat.getPairwise(nodeI.posNum, rotI.labelNum, nodeJ.posNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotJE = -this.emat.getOneBody(nodeJ.posNum, rotJ.labelNum);
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
                        }
                    }
                    //Now we send message from i->j
                    double partFunctionJ = 0.0;
                    for (int rotJindex = 0; rotJindex < nodeJ.labelList.size(); rotJindex++) {
                        double messageAtRotJ = 0.0;
                        for (int rotIindex = 0; rotIindex < nodeI.labelList.size(); rotIindex++) {
                            MRFLabel rotI = nodeI.labelList.get(rotIindex);
                            MRFLabel rotJ = nodeJ.labelList.get(rotJindex);

                            double pairwiseE = -this.emat.getPairwise(nodeI.posNum, rotI.labelNum, nodeJ.posNum, rotJ.labelNum);
                            double edgeProbability = getEdgeProbability(i, j);
                            double rotIE = -this.emat.getOneBody(nodeI.posNum, rotI.labelNum);
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
                        }
                    }
                }
            }
        }
        this.maxChange = maxChangeMessage;
        return updatedMessages;
    }

    private void checkMessages(double[][][] messages) {
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                double sumJ = 0.0;
                for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                    sumJ += messages[nodeI.posNum][nodeJ.posNum][rotJ];
                }
                if (Math.abs(1 - sumJ) > 1e-6) {
                    throw new RuntimeException("Messages Not Normalized between node " + nodeI.posNum + " and node " + nodeJ.posNum);
                }
                double sumI = 0.0;
                for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                    sumI += messages[nodeJ.posNum][nodeI.posNum][rotI];
                }
                if (Math.abs(1 - sumI) > 1e-6) {
                    throw new RuntimeException("Messages Not Normalized between node " + nodeJ.posNum + " and node " + nodeI.posNum);
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
                MRFNode nodeI = this.nodeList.get(i);
                MRFNode nodeJ = this.nodeList.get(j);
                updateEdgeMarginal(nodeI, nodeJ, messages);
            }
        }
    }

    private void updateEdgeMarginal(MRFNode nodeI, MRFNode nodeJ, double[][][] messages) {
        double partFuncOverall = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double pairwiseE = this.emat.getPairwise(nodeI.posNum, nodeI.labelList.get(rotI).labelNum, nodeJ.posNum, nodeJ.labelList.get(rotJ).labelNum);
                double edgeProb = getEdgeProbability(nodeI.posNum, nodeJ.posNum);
                double nodeIE = this.emat.getOneBody(nodeI.posNum, nodeI.labelList.get(rotI).labelNum);
                double nodeJE = this.emat.getOneBody(nodeJ.posNum, nodeJ.labelList.get(rotJ).labelNum);

//                double normalized = (-((pairwiseE / edgeProb) + nodeIE + nodeJE)); //- expNorm;
                double normalized = (-((pairwiseE / edgeProb) + nodeIE + nodeJE)) - this.expNormMarginals[nodeI.posNum][nodeJ.posNum];
                double marginal = Math.exp(normalized / this.constRT);
//                double marginal = Math.exp(normalized);
                marginal = marginal * (getProductMessages(nodeI, rotI, nodeJ, messages));
                marginal = marginal * (getProductMessages(nodeJ, rotJ, nodeI, messages));

                partFuncOverall += marginal;
                this.marginalProbabilies.setPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ, marginal);
            }
        }
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double unNormalized = this.marginalProbabilies.getPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ);
                double normalized = unNormalized / partFuncOverall;
                this.marginalProbabilies.setPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ, normalized);
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
            sum += this.marginalProbabilies.getOneBody(node.posNum, rot);
        }
        if (Math.abs(sum - 1) > 1e-6) {
            throw new RuntimeException("Marginal Not Normalized at Node: " + node.posNum);
        }
    }

    private void checkPairwiseMarginal(MRFNode nodeI, MRFNode nodeJ) {
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            double marginal = this.marginalProbabilies.getOneBody(nodeI.posNum, rotI);
            double sum = 0.0;
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                sum += this.marginalProbabilies.getPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ);
            }
            if (Math.abs(sum - marginal) > 0.0001) {
                throw new RuntimeException("Marginal Not Normalized at Node Pairs: " + nodeI.posNum + " " + nodeJ.posNum);
            }
        }
        for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
            double marginal = this.marginalProbabilies.getOneBody(nodeJ.posNum, rotJ);
            double sum = 0.0;
            for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                sum += this.marginalProbabilies.getPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ);
            }
            if (Math.abs(sum - marginal) > 0.0001) {
                throw new RuntimeException("Marginal Not Normalized at Node Pairs: " + nodeJ.posNum + " " + nodeI.posNum);
            }
        }
    }

    private void updateNodeMarginal(MRFNode node, double[][][] messages) {
        double partFunc = 0.0;
        double expNorm = Double.NEGATIVE_INFINITY;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            expNorm = Math.max(expNorm, -this.emat.getOneBody(node.posNum, node.labelList.get(rot).labelNum));
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double normalized = -this.emat.getOneBody(node.posNum, node.labelList.get(rot).labelNum) - expNorm;
            double update = Math.exp(normalized / this.constRT);
//            double update = Math.exp(normalized);
            for (MRFNode neighbor : node.neighborList) {
                double message = messages[neighbor.posNum][node.posNum][rot];
                double prob = getEdgeProbability(neighbor.posNum, node.posNum);
                update = update * (Math.pow(message, prob));
            }
            this.marginalProbabilies.setOneBody(node.posNum, rot, update);
            partFunc += update;
        }
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double unNormalized = this.marginalProbabilies.getOneBody(node.posNum, rot);
            this.marginalProbabilies.setOneBody(node.posNum, rot, unNormalized / partFunc);
        }
    }

    private double getProductMessages(MRFNode nodeReceiving, int labelIndex, MRFNode nodeExcluded, double[][][] currentMessages) {
        double product = 1;
        for (MRFNode node : nodeReceiving.neighborList) {
            if (!(node.equals(nodeExcluded))) {
                double message = currentMessages[node.posNum][nodeReceiving.posNum][labelIndex];
                double probability = getEdgeProbability(nodeReceiving.posNum, node.posNum);
                product = product * (Math.pow(message, probability));
            }
        }
        double messageNodeExcluded = currentMessages[nodeExcluded.posNum][nodeReceiving.posNum][labelIndex];
        double probability = getEdgeProbability(nodeExcluded.posNum, nodeReceiving.posNum);
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
                    messages[nodeI.posNum][nodeJ.posNum][rotJ] = 1.0 / (double) nodeJ.labelList.size();
                }
                for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
                    messages[nodeJ.posNum][nodeI.posNum][rotI] = 1.0 / (double) nodeI.labelList.size();
                }
                //Then normalize j-i x_i
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
            double E = this.emat.getOneBody(node.posNum, node.labelList.get(rot).labelNum);
            double prob = this.marginalProbabilies.getOneBody(node.posNum, rot);
            enthalpy += E * prob;
        }
        return enthalpy;
    }

    private double getPairwiseNodeEnthalpy(MRFNode nodeI, MRFNode nodeJ) {
        double enthalpy = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double E = emat.getPairwise(nodeI.posNum, nodeI.labelList.get(rotI).labelNum, nodeJ.posNum, nodeJ.labelList.get(rotJ).labelNum);
                double prob = this.marginalProbabilies.getPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ);
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
                if (this.interactionGraph[node1.posNum][node2.posNum]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                }
            }
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (int rot = 0; rot < node.labelList.size(); rot++) {
            double prob = this.marginalProbabilies.getOneBody(node.posNum, rot);
            entropy += (-1.0) * prob * Math.log(prob);
        }
        return entropy;
    }

    private double getMutualInformation(MRFNode nodeI, MRFNode nodeJ) {
        double mutualInf = 0.0;
        for (int rotI = 0; rotI < nodeI.labelList.size(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.labelList.size(); rotJ++) {
                double probIJ = this.marginalProbabilies.getPairwise(nodeI.posNum, rotI, nodeJ.posNum, rotJ);
                double probI = this.marginalProbabilies.getOneBody(nodeI.posNum, rotI);
                double probJ = this.marginalProbabilies.getOneBody(nodeJ.posNum, rotJ);
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

    double[][] getEdgeWeights() {
        double[][] edgeWeights = new double[this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            double[] edgeFromI = new double[i];
            MRFNode nodeI = this.nodeList.get(i);
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (interactionGraph[i][j]) {
                    edgeFromI[j] = -getMutualInformation(nodeI, nodeJ);
                }
            }
            edgeWeights[i] = edgeFromI;
        }
        return edgeWeights;
    }

    private void updateEdgeProbabilies(double[][] descentDirection, int iteration) {
        double stepSize = 2.0 / (iteration + 4.0);
        for (int i = 0; i < this.numNodes; i++) {
            for (int j = 0; j < i; j++) {
                this.edgeProbabilities[i][j] = stepSize * descentDirection[i][j] + (1 - stepSize) * edgeProbabilities[i][j];
            }
        }
    }

    public double getLogZ() {
        return this.logZ;
    }
}
