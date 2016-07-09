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
public class TRBP {

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
//    static int maxNumEdgeUpdates = 3;
    static int maxNumEdgeUpdates = 3;

    public double[][] edgeProbabilities;
    //Edge Weights are given by the negative mutual information
    //They are used for to update edge probabilities via the MST
    double[][] edgeWeights;

    //The log messages
    double[][][] logMessages;

    // We precompute everything in the exponential on equation 39.
    // These are updated everytime we update edge probabilities
    double[][][][] messageEnergies;

    // This is 2*number of edges
    int numMessagesPairs;

    public ExpFunction ef = new ExpFunction();

    double averageChange;

    //This records the accuracy with which we want to continue updating messages
    double accuracyWithinEdgeProb = 0.001;
    double messageConvergence = 1e-5;
    double accuracyBetweenEdgeProb = 0.001;

    // When used in the branch and bound algorithm, we can save time not updating 
    // edges if the child partition function is already much smaller than the parent
    double parentUpperBound = Double.NEGATIVE_INFINITY;
    boolean useParentUpperBound = false;
    double cutOffUpperBound = 3; // Cutoff which we can skip edge updates

    static public boolean verbose = true;

    // Node weights are used in the branch and bound algorithm
    public double[] nodeWeights;

    public TRBP(MarkovRandomField mrf) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.interactionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.getNumLabels();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.interactionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);
        this.numMessagesPairs = 2 * getNumEdges(this.interactionGraph);

        initializeEdgeWeights();

        runTRBP();
    }

    public TRBP(MarkovRandomField mrf, double parentUpperBound) {
        this.nodeList = mrf.nodeList;
        this.emat = mrf.emat;
        this.interactionGraph = mrf.interactionGraph;
        this.numNodes = nodeList.size();

        this.numLabelsPerNode = new int[this.numNodes];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode node = this.nodeList.get(i);
            numLabelsPerNode[i] = node.getNumLabels();
        }

        this.edgeProbabilities = initializeEdgeProbabilities(this.interactionGraph);
        this.logMessages = initializeLogMessages(0.0);

        this.marginalProbabilities = new TupleMatrix(numNodes, numLabelsPerNode, Double.POSITIVE_INFINITY, 0.0);
        this.numMessagesPairs = 2 * getNumEdges(this.interactionGraph);

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
            this.messageEnergies = precomputeMessageEnergies();

            //Keep track of last logZ from previous message update
            double lastLogZMessage = Double.POSITIVE_INFINITY;
            int numMessageUpdates = 0;
            double changeBetweenMessageUpdates = Double.POSITIVE_INFINITY;
            // While we have not converged...
            while (((this.averageChange > this.messageConvergence) || changeBetweenMessageUpdates > this.accuracyWithinEdgeProb) || numMessageUpdates < 10) {
                this.averageChange = 0.0;
                //This is the meat of the algorithm
                updateMessagesParallel(logMessages, damping);
//                updateMessagesSequentially(logMessages, useDamping);
                // Update the marginals and compute LogZ every 10 runs
                if (numMessageUpdates % 10 == 0 && (this.averageChange < this.messageConvergence)) {
                    updateMarginals();
                    double currentlogZ = calcUBLogZ();
                    // Get the change in logZ
                    // Classically, people use change in messages rather than change
                    // in logZ, which may be better... we use both
                    double change = Math.abs(lastLogZMessage - currentlogZ);
                    lastLogZMessage = currentlogZ;
                    changeBetweenMessageUpdates = change;
                    if (verbose) {
                        System.out.println("Average Message Change: " + this.averageChange);
                    }
                } else {
                    if (verbose) {
                        System.out.println("Average Message Change: " + this.averageChange);
                        /*                        updateMarginals();
                         double currentLogZ = calcUBLogZ();
                         System.out.println("UBLogZ: " + currentLogZ);*/
                    }
                }
                numMessageUpdates++;
            }
//            updateMarginals();
//            double currentlogZ = calcUBLogZ();
            lastLogZEdge = lastLogZMessage;
            this.logZ = Math.min(this.logZ, lastLogZEdge);
            // This is (only) useful for the partFuncTree calculation
            // Basically, we don't need high accuracy if the partition function is very small
            // Compared to the parent node 
            if (useParentUpperBound && (lastLogZEdge + cutOffUpperBound < parentUpperBound)) {
                break;
            }
            numEdgeUpdates++;
        }
        double[] degrees = GraphUtils.getWeightedDegrees(edgeWeights);
        if (degrees == null) {
            throw new RuntimeException("Degrees is Null");
        }
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

    double[][][][] precomputeMessageEnergies() {
        double[][][][] messageTS = new double[this.numNodes][this.numNodes][][];
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode nodeI = this.nodeList.get(i);
            for (int j = 0; j < this.numNodes; j++) {
                if (this.interactionGraph[i][j]) {
                    MRFNode nodeJ = this.nodeList.get(j);
                    double edgeProb = getEdgeProbability(nodeI, nodeJ);

                    messageTS[i][j] = new double[nodeJ.getNumLabels()][nodeI.getNumLabels()];
                    for (int labelJNum = 0; labelJNum < nodeJ.getNumLabels(); labelJNum++) {
                        int labelJ = nodeJ.labels[labelJNum];
                        for (int labelINum = 0; labelINum < nodeI.getNumLabels(); labelINum++) {
                            int labelI = nodeI.labels[labelINum];

                            double pairPot = getPairwisePotential(nodeI, labelI, nodeJ, labelJ);
                            double singlePotI = getOneBodyPotential(nodeI, labelI);
                            messageTS[i][j][labelJNum][labelINum] = ((pairPot / edgeProb) + singlePotI) / this.constRT;
                        }
                    }
                }
            }
        }
        return messageTS;
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

        int[][] messageOrdering = getParallelMessagePassingOrdering();
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
                sumLogMessages = getSumLogMessages(nodeI);
                previousSendingNodeIndex = nodeI.index;
            }

            //Get updated log Messages from nodeI to nodeJ
            double[] logMessagesNodeIToNodeJ = getUpdatedLogMessage(nodeI, nodeJ, sumLogMessages, this.messageEnergies[nodePair[0]][nodePair[1]]);

            double[] previousLogMessages = previousMessages[nodeI.index][nodeJ.index];

            double change = dampenMessages(logMessagesNodeIToNodeJ, previousLogMessages, nodeI.index, nodeJ.index, updatedLogMessages, damping);
            this.averageChange += change;
        }

        //Update average change
        this.averageChange = this.averageChange / this.numMessagesPairs;
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
            double change;
            if (Double.isFinite(damped) && Double.isFinite(prevMessage)) {
                change = Math.abs(damped - prevMessage);
            }
            else if (Double.isInfinite(damped) && Double.isInfinite(prevMessage)){
                change = 0.0;
            }
            else {
                change = Double.POSITIVE_INFINITY;
            }
            sumChange += change;
        }
        return sumChange / (double) updated.length;
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
        double[] marginals = new double[node.getNumLabels()];
        for (int labelIndex = 0; labelIndex < node.getNumLabels(); labelIndex++) {
            int label = node.labels[labelIndex];

            double labelPot = getOneBodyPotential(node, label) / this.constRT;

            double sumWeightedLogMessages = 0.0;
            for (MRFNode neighbor : this.nodeList) {
                if (this.interactionGraph[node.index][neighbor.index]) {
                    double edgeProb = getEdgeProbability(node, neighbor);
                    double logMessageNeighborToNode = getLogMessage(neighbor, node, labelIndex);
                    sumWeightedLogMessages += edgeProb * logMessageNeighborToNode;
                }
            }
            if (Double.isNaN(sumWeightedLogMessages)){
                throw new RuntimeException("sumWeightedLogMessages is NaN");
            }
            normalizer = Math.max(normalizer, labelPot + sumWeightedLogMessages);
            marginals[labelIndex] = labelPot + sumWeightedLogMessages;
            if (Double.isNaN(marginals[labelIndex])){
                throw new RuntimeException("Marginal is NaN");
            }
        }
        for (int labelIndex = 0; labelIndex < node.getNumLabels(); labelIndex++) {
            marginals[labelIndex] -= normalizer;
            marginals[labelIndex] = Math.exp(marginals[labelIndex]);
            partFunc += marginals[labelIndex];
        }
        for (int labelIndex = 0; labelIndex < node.getNumLabels(); labelIndex++) {
            double nonNormalizedMarginal = marginals[labelIndex];
            this.marginalProbabilities.setOneBody(node.index, labelIndex, nonNormalizedMarginal / partFunc);
            if (Double.isNaN(this.marginalProbabilities.getOneBody(node.index, labelIndex))){
                System.out.println("nonNormMarg: "+nonNormalizedMarginal);
                System.out.println("partFunc: "+partFunc);
                throw new RuntimeException("marginal is NaN");
            }
        }
    }

    private void updateEdgeMarginal(MRFNode nodeI, MRFNode nodeJ) {
        double partFunc = 0;
        double normalizer = Double.NEGATIVE_INFINITY;
        double[] sumLogMessagesItoJ = getSumLogMessages(nodeI);
        sumLogMessagesItoJ = subtractReceiver(sumLogMessagesItoJ, nodeI, nodeJ);

        double[] sumLogMessagesJtoI = getSumLogMessages(nodeJ);
        sumLogMessagesJtoI = subtractReceiver(sumLogMessagesJtoI, nodeJ, nodeI);

        for (int labelIndexI = 0; labelIndexI < nodeI.getNumLabels(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.getNumLabels(); labelIndexJ++) {
                int labelI = nodeI.labels[labelIndexI];
                int labelJ = nodeJ.labels[labelIndexJ];
                double pairPot = getPairwisePotential(nodeI, labelI, nodeJ, labelJ);
                double edgeProb = getEdgeProbability(nodeI, nodeJ);
                double labelIPot = getOneBodyPotential(nodeI, labelI);
                double labelJPot = getOneBodyPotential(nodeJ, labelJ);

                double potential = ((pairPot / edgeProb) + labelIPot + labelJPot) / this.constRT;
                double sumWeightedLogMessages = sumLogMessagesItoJ[labelIndexI] + sumLogMessagesJtoI[labelIndexJ];
                normalizer = Math.max(normalizer, potential + sumWeightedLogMessages);

                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, potential + sumWeightedLogMessages);
            }
        }
        for (int labelIndexI = 0; labelIndexI < nodeI.getNumLabels(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.getNumLabels(); labelIndexJ++) {
                double unNormalized = this.marginalProbabilities.getPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ);
                double normalized = unNormalized - normalizer;
                double exponentiated = Math.exp(normalized);
                this.marginalProbabilities.setPairwise(nodeI.index, labelIndexI, nodeJ.index, labelIndexJ, exponentiated);
                partFunc += exponentiated;
            }
        }
        for (int labelIndexI = 0; labelIndexI < nodeI.getNumLabels(); labelIndexI++) {
            for (int labelIndexJ = 0; labelIndexJ < nodeJ.getNumLabels(); labelIndexJ++) {
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
    int[][] getSequentialMessagePassingOrdering(boolean sequential) {
        int[][] messagePassingOrdering = new int[numMessagesPairs][];

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
    int[][] getParallelMessagePassingOrdering() {
        int[][] messagePassingOrdering = new int[numMessagesPairs][];

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

    double[] getUpdatedLogMessage(MRFNode sendingNode, MRFNode receivingNode, double[] sumLogMessages, double[][] potentials) {
        double[] updatedLogMessages = new double[receivingNode.getNumLabels()];
        double largestLogMessage = Double.NEGATIVE_INFINITY;

//        ArrayList<Double> oneBodyESending = this.emat.getOneBodyEnergies(sendingNode.posNum);
//        ArrayList<ArrayList<Double>> twoBodyE = this.emat.getTwoBodyEnergies(sendingNode.posNum, receivingNode.posNum);

        /*double[] sumLogMessagesMinusReceiver = new double[sendingNode.getNumLabels()];
         for (int i = 0; i < sumLogMessages.length; i++) {
         sumLogMessagesMinusReceiver[i] = sumLogMessages[i] - getLogMessage(receivingNode, sendingNode, sendingNode.labelList.get(i));
         }*/
        double[] sumLogMessagesMinusReceiver = subtractReceiver(sumLogMessages, sendingNode, receivingNode);
        for (int receivingLabelIndex = 0; receivingLabelIndex < receivingNode.getNumLabels(); receivingLabelIndex++) {
//            int receivingLabel = receivingNode.labels[receivingLabelIndex];

            double normalizer = Double.NEGATIVE_INFINITY;
            double[] toBeSummed = new double[sendingNode.getNumLabels()];

            double[] potentialAtReceiveLabel = potentials[receivingLabelIndex];
            for (int sendingLabelIndex = 0; sendingLabelIndex < sendingNode.getNumLabels(); sendingLabelIndex++) {
//                int sendingLabel = sendingNode.labels[sendingLabelIndex];

                //compute everything between the curly brackets in eq 39. of TRBP Paper
                //We do this in the log domain and exponentiate aftern
//                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel);
                // Optimization: we precompute these potentials (everything in the exponential in eq. 39)
                double pot = potentialAtReceiveLabel[sendingLabelIndex];
//                double pairwisePot = getPairwisePotential(sendingNode, sendingLabel, receivingNode, receivingLabel, twoBodyE);
//                double edgeProbability = getEdgeProbability(receivingNode, sendingNode);
//                double sendingLabelPot = getOneBodyPotential(sendingNode, sendingLabel);
//                double sendingLabelPot = getOneBodyPotential(sendingLabel, oneBodyESending);

//                double sumLogMessage = sumLogMessages[sendingLabelIndex] - getLogMessage(receivingNode, sendingNode, sendingLabel);
                double sumLogMessage = sumLogMessagesMinusReceiver[sendingLabelIndex];
                if (Double.isNaN(sumLogMessage)) {
                    throw new RuntimeException("SumLogMessage is Nan");
                }
//                toBeSummed[sendingLabelIndex] = (((pairwisePot / edgeProbability) + sendingLabelPot) / constRT) + sumLogMessage;
                //toBeSummed[sendingLabelIndex] = getMessageUpdateInSum(pairwisePot, edgeProbability, sendingLabelPot, sumLogMessage);
                toBeSummed[sendingLabelIndex] = getMessageUpdateInSum(pot, sumLogMessage);

                if (Double.isNaN(toBeSummed[sendingLabelIndex])) {
                    throw new RuntimeException("TOBESUMMED NAN");
                }

                normalizer = Math.max(normalizer, toBeSummed[sendingLabelIndex]);
            }

            int numToBeSummed = sendingNode.getNumLabels();
            double sum;
            if (Double.isInfinite(normalizer)) {
                sum = sumNormalizedToBeSummed(toBeSummed, 0.0, numToBeSummed);
            } else {
                sum = sumNormalizedToBeSummed(toBeSummed, normalizer, numToBeSummed);
            }
            if (Double.isNaN(sum)) {
                throw new RuntimeException("Sum is Nan");
            }
//            double nonNormalizedLogMessage = Math.log(sum) + normalizer;
            double nonNormalizedLogMessage = getNonNormalizedLogMessage(sum, normalizer);

            if (Double.isNaN(nonNormalizedLogMessage)) {
                throw new RuntimeException("nonNormlizedLogM is Nan");
            }
            largestLogMessage = Math.max(largestLogMessage, nonNormalizedLogMessage);
            //int messageIndex = receivingNode.labelList.indexOf(receivingLabel);
            updatedLogMessages[receivingLabelIndex] = nonNormalizedLogMessage;
        }
        if (Double.isInfinite(largestLogMessage)) {
            return updatedLogMessages;
            //   throw new RuntimeException("LARGEST LOG MESSAGE INFINITE");
        }
        
        //Normalize by subtracting the largest logMessage;
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
            double logMessageToSubtract = getLogMessage(receiver, sender, i);
            if (Double.isFinite(logMessageToSubtract)) {
                sumLogMessagesMinusReceiver[i] = sumLogMessages[i] - getLogMessage(receiver, sender, i);
            } else {
//                System.out.println("Using new technique");
                double sum = 0.0;
                for (int j = 0; j < sender.getNumLabels(); j++) {
                    if (j != i) {
                        sum += getLogMessage(receiver, sender, j);
                    }
                }
                sumLogMessagesMinusReceiver[i] = sum;
            }
//            System.out.println("sumLogMessages "+sumLogMessages[i]);
//            System.out.println("logMessageMinus "+getLogMessage(receiver, sender, i));
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

    public double getMessageUpdateInSum(double pot, double sumLogMessage) {
        return pot + sumLogMessage;
    }

    double getSumLogMessage(MRFNode receivingNode, int labelIndex) {
        double sum = 0;
        int receivingIndex = receivingNode.index;
        for (int neighborIndex = 0; neighborIndex < this.numNodes; neighborIndex++) {
            if (this.interactionGraph[receivingIndex][neighborIndex]) {
                MRFNode nodeV = this.nodeList.get(neighborIndex);
                double edgeProbVToSender = getEdgeProbability(receivingNode, nodeV);
                double logMessageVToSender = getLogMessage(nodeV, receivingNode, labelIndex);
                sum += edgeProbVToSender * logMessageVToSender;
            }
        }
        return sum;
    }

    double[] getSumLogMessages(MRFNode sendingNode) {
        double[] sumLogMessages = new double[sendingNode.getNumLabels()];
        int index = 0;
        for (int labelIndex = 0; labelIndex < sendingNode.getNumLabels(); labelIndex++) {
            sumLogMessages[index] = getSumLogMessage(sendingNode, labelIndex);
            index++;
        }
        return sumLogMessages;
    }

    public int getLabelListSize(MRFNode node) {
        return node.getNumLabels();
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
            if (Double.isNaN(enthalpy)){
                throw new RuntimeException("Enthalpy is NaN after Single Node");
            }
            for (int j = 0; j < i; j++) {
                MRFNode node2 = nodeList.get(j);
                if (this.interactionGraph[node1.index][node2.index]) {
                    enthalpy += getPairwiseNodeEnthalpy(node1, node2);
                    if (Double.isNaN(enthalpy)){
                        throw new RuntimeException("Enthalpy is NaN after Double Node");
                    }
                }
            }
        }
        if (Double.isNaN(enthalpy)){
            throw new RuntimeException("Enthalpy is NaN");
        }
        return enthalpy;
    }

    public double getEntropy() {
        double entropy = 0.0;
        for (int i = 0; i < this.numNodes; i++) {
            MRFNode nodeI = this.nodeList.get(i);
            double singleNodeEntropy = getSingleNodeEntropy(nodeI);
            entropy += singleNodeEntropy;
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                if (interactionGraph[i][j]) {
                    double edgeProb = getEdgeProbability(nodeI, nodeJ);
                    double mutualInf = getMutualInformation(nodeI, nodeJ);
                    //Update Edge Weights to Save Computational Cost
                    this.edgeWeights[i][j] = -mutualInf;///Math.sqrt(entropyI*entropyJ);
                    entropy -= edgeProb * mutualInf;
                }
            }
        }
        if (Double.isNaN(entropy)){
            throw new RuntimeException("Entropy is NaN");
        }
        return entropy;
    }

    private double getSingleNodeEnthalpy(MRFNode node) {
        double enthalpy = 0.0;
        for (int i = 0; i < node.getNumLabels(); i++) {
            int label = node.labels[i];
            double E = this.emat.getOneBody(node.posNum, label);
            double prob = this.marginalProbabilities.getOneBody(node.index, i);
            if (Double.isNaN(E)){
                throw new RuntimeException("E is NaN in SNEnth");
            }
            if (Double.isNaN(prob)){
                throw new RuntimeException("prob is NaN in SNEnth");
            }
            enthalpy += E * prob;
        }
        return enthalpy;
    }

    private double getSingleNodeEntropy(MRFNode node) {
        double entropy = 0.0;
        for (int rot = 0; rot < node.getNumLabels(); rot++) {
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
        for (int rotI = 0; rotI < nodeI.getNumLabels(); rotI++) {
            for (int rotJ = 0; rotJ < nodeJ.getNumLabels(); rotJ++) {
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
        for (int i = 0; i < nodeI.getNumLabels(); i++) {
            int labelI = nodeI.labels[i];
            for (int j = 0; j < nodeJ.getNumLabels(); j++) {
                int labelJ = nodeJ.labels[j];
                double E = emat.getPairwise(nodeI.posNum, labelI, nodeJ.posNum, labelJ);
                double prob = this.marginalProbabilities.getPairwise(nodeI.index, i, nodeJ.index, j);
                if (Double.isInfinite(E)) {
                    enthalpy += 0.0;
                } else {
                    enthalpy += E * prob;
                }
            }
        }
        return enthalpy;
    }

    private void initializeEdgeWeights() {
        this.edgeWeights = new double[this.numNodes][];
        for (int i = 0; i < this.numNodes; i++) {
            this.edgeWeights[i] = new double[i];
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

    double getLogMessage(MRFNode sendingNode, MRFNode receivingNode, int receivingLabelIndex
    ) {
        return this.logMessages[sendingNode.index][receivingNode.index][receivingLabelIndex];
    }

    double getPairwisePotential(MRFNode nodeI, int labelI, MRFNode nodeJ, int labelJ
    ) {
        return -this.emat.getPairwise(nodeI.posNum, labelI, nodeJ.posNum, labelJ);
    }

    double getPairwisePotential(MRFNode nodeI, int labelI, MRFNode nodeJ, int labelJ, ArrayList<ArrayList<Double>> energies) {
        if (nodeI.posNum > nodeJ.posNum) {
            return -energies.get(labelI).get(labelJ);
        }
        return -energies.get(labelJ).get(labelI);
    }

    double getOneBodyPotential(MRFNode node, int label
    ) {
        return -this.emat.getOneBody(node.posNum, label);
    }

    double getOneBodyPotential(int label, ArrayList<Double> energies) {
        return -energies.get(label);
    }

    public double getLogZ() {
        return this.logZ;
    }

    public static void setNumEdgeProbUpdates(int numUpdates) {
        assert (numUpdates >= 0);
        maxNumEdgeUpdates = numUpdates;
    }

}
