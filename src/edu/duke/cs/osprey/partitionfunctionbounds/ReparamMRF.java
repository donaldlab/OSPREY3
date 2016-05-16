/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class ReparamMRF {

    public ArrayList<MRFNode> nodeList = new ArrayList<>();
    ArrayList<MRFNode> clampedNodeList = new ArrayList<>();
    boolean[][] interactionGraph;
    boolean[][] nonClampedInteractionGraph;

    int numPos;

    UpdatedEmat emat;
    PruningMatrix pruneMat;

    public ReparamMRF(SearchProblem searchProblem, double eCut) {
        this.pruneMat = searchProblem.pruneMat;

        ConfSpace cSpace = searchProblem.confSpace;
        this.numPos = cSpace.numPos;

        //create nodeList
        int numClamped = 0;
        int numNonClamped = 0;
        ArrayList<MRFNode> allNodes = new ArrayList();
        for (int pos = 0; pos < numPos; pos++) {
            ArrayList<Integer> unprunedRCs = pruneMat.unprunedRCsAtPos(pos);
            if (unprunedRCs.size() > 1) {
                MRFNode node = new MRFNode(pos, unprunedRCs, numNonClamped);
                nodeList.add(node);
                numNonClamped++;

                allNodes.add(node);
            } else {
                //Node is clamped
                MRFNode node = new MRFNode(pos, unprunedRCs, numClamped);
                clampedNodeList.add(node);
                numClamped++;

                allNodes.add(node);
            }
        }

        //create interaction graph
        createEnergyInteractionGraph(searchProblem.emat, eCut, allNodes);
        createNonClampedInteractionGraph();
        //create neighborList for each node
        for (MRFNode node : this.nodeList) {
            node.neighborList = getNeighbors(node, this.interactionGraph);
        }

        this.emat = new UpdatedEmat(searchProblem.emat, clampedNodeList, interactionGraph);
    }

    /**
     * Returns an array that maps the index into the (non-clamped) nodeList to
     * the pos num. This is useful in the branch-and-bound algorithms
     *
     * @return
     */
    public int[] getIndexToPosNumMap() {
        int[] indexToPosNum = new int[this.nodeList.size()];
        for (int index = 0; index < this.nodeList.size(); index++) {
            int posNum = this.nodeList.get(index).posNum;
            indexToPosNum[index] = posNum;
        }
        return indexToPosNum;
    }

    public int getNumEdge() {
        int numEdges = 0;
        for (int i = 0; i < this.nodeList.size(); i++) {
            for (int j = 0; j < i; j++) {
                if (this.nonClampedInteractionGraph[i][j]){
                    numEdges++;
                }
            }
        }
        return numEdges;
    }

    public ReparamMRF(EnergyMatrix aemat, PruningMatrix apruneMat, double eCut) {
        this.pruneMat = apruneMat;

        this.numPos = aemat.numPos();

        //create nodeList
        int numClamped = 0;
        int numNonClamped = 0;
        ArrayList<MRFNode> allNodes = new ArrayList();
        for (int pos = 0; pos < numPos; pos++) {
            ArrayList<Integer> unprunedRCs = pruneMat.unprunedRCsAtPos(pos);
            if (unprunedRCs.size() > 1) {
                MRFNode node = new MRFNode(pos, unprunedRCs, numNonClamped);
                nodeList.add(node);
                numNonClamped++;

                allNodes.add(node);
            } else {
                //Node is clamped
                MRFNode node = new MRFNode(pos, unprunedRCs, numClamped);
                clampedNodeList.add(node);
                numClamped++;

                allNodes.add(node);
            }
        }

        //create interaction graph
        createEnergyInteractionGraph(aemat, eCut, allNodes);
        createNonClampedInteractionGraph();
        //create neighborList for each node
        for (MRFNode node : this.nodeList) {
            node.neighborList = getNeighbors(node, this.interactionGraph);
        }
        this.emat = new UpdatedEmat(aemat, clampedNodeList, interactionGraph);
    }

    public ReparamMRF(SearchProblem searchProblem, int[] partialNode, double eCut) {
        this.pruneMat = searchProblem.pruneMat;

        ConfSpace cSpace = searchProblem.confSpace;
        this.numPos = cSpace.numPos;

        //create nodeList
        int numClamped = 0;
        int numNonClamped = 0;
        ArrayList<MRFNode> allNodes = new ArrayList();
        for (int pos = 0; pos < numPos; pos++) {
            ArrayList<Integer> unprunedRCs;
            if (partialNode[pos] == -1) {
                unprunedRCs = pruneMat.unprunedRCsAtPos(pos);
            } else {
                unprunedRCs = new ArrayList();
                unprunedRCs.add(partialNode[pos]);
            }

            if (unprunedRCs.size() > 1) {
                MRFNode node = new MRFNode(pos, unprunedRCs, numNonClamped);
                nodeList.add(node);
                numNonClamped++;

                allNodes.add(node);
            } else {
                //Node is clamped
                MRFNode node = new MRFNode(pos, unprunedRCs, numClamped);
                clampedNodeList.add(node);
                numClamped++;

                allNodes.add(node);
            }
        }

        //create interaction graph
        createEnergyInteractionGraph(searchProblem.emat, eCut, allNodes);
        createNonClampedInteractionGraph();
        //create neighborList for each node
        for (MRFNode node : this.nodeList) {
            node.neighborList = getNeighbors(node, this.interactionGraph);
        }

        this.emat = new UpdatedEmat(searchProblem.emat, clampedNodeList, interactionGraph);
    }

    public ReparamMRF(EnergyMatrix aemat, PruningMatrix apruneMat, int[] partialNode, double eCut) {
        this.pruneMat = apruneMat;

        this.numPos = aemat.numPos();

        //create nodeList
        int numClamped = 0;
        int numNonClamped = 0;
        ArrayList<MRFNode> allNodes = new ArrayList();
        for (int pos = 0; pos < numPos; pos++) {
            ArrayList<Integer> unprunedRCs;
            if (partialNode[pos] == -1) {
                unprunedRCs = pruneMat.unprunedRCsAtPos(pos);
            } else {
                unprunedRCs = new ArrayList();
                unprunedRCs.add(partialNode[pos]);
            }

            if (unprunedRCs.size() > 1) {
                MRFNode node = new MRFNode(pos, unprunedRCs, numNonClamped);
                nodeList.add(node);
                numNonClamped++;

                allNodes.add(node);
            } else {
                //Node is clamped
                MRFNode node = new MRFNode(pos, unprunedRCs, numClamped);
                clampedNodeList.add(node);
                numClamped++;

                allNodes.add(node);
            }
        }

        //create interaction graph
        createEnergyInteractionGraph(aemat, eCut, allNodes);
        createNonClampedInteractionGraph();
        //create neighborList for each node
        for (MRFNode node : this.nodeList) {
            node.neighborList = getNeighbors(node, this.interactionGraph);
        }

        this.emat = new UpdatedEmat(aemat, clampedNodeList, interactionGraph);
    }

    private void createNonClampedInteractionGraph() {
        int numNodes = this.nodeList.size();
        this.nonClampedInteractionGraph = new boolean[numNodes][numNodes];
        for (int i = 0; i < numNodes; i++) {
            MRFNode nodeI = this.nodeList.get(i);
            nonClampedInteractionGraph[i][i] = false;
            for (int j = 0; j < i; j++) {
                MRFNode nodeJ = this.nodeList.get(j);
                nonClampedInteractionGraph[i][j] = this.interactionGraph[nodeI.posNum][nodeJ.posNum];
                nonClampedInteractionGraph[j][i] = this.interactionGraph[nodeJ.posNum][nodeI.posNum];
            }
        }
    }

    private void createEnergyInteractionGraph(EnergyMatrix aemat, double eCut, ArrayList<MRFNode> allNodes) {
        this.interactionGraph = new boolean[numPos][numPos];
        int countInteraction = 0;
        int possibleInteraction = 0;

        if (eCut == 0.0) {
            //initialize all values to true
            for (int nodeNum1 = 0; nodeNum1 < numPos; nodeNum1++) {
                for (int nodeNum2 = 0; nodeNum2 < numPos; nodeNum2++) {
                    if (nodeNum1 == nodeNum2) {
                        interactionGraph[nodeNum1][nodeNum2] = false;
                    } else {
                        interactionGraph[nodeNum1][nodeNum2] = true;
                    }
                }
            }
        }

        //otherwise initialize all values to false
        for (int nodeNum1 = 0; nodeNum1 < numPos; nodeNum1++) {
            for (int nodeNum2 = 0; nodeNum2 < numPos; nodeNum2++) {
                interactionGraph[nodeNum1][nodeNum2] = false;
            }
        }
        //now get maxInteraction and check if it is greater than eCut
        for (int nodeNum1 = 0; nodeNum1 < numPos; nodeNum1++) {
            for (int nodeNum2 = 0; nodeNum2 < nodeNum1; nodeNum2++) {
                double maxInteraction = 0.0;
                MRFNode node1 = allNodes.get(nodeNum1);
                MRFNode node2 = allNodes.get(nodeNum2);
                for (int label1 : node1.labels) {
                    for (int label2 : node2.labels) {
                        double pairE = aemat.getPairwise(node1.posNum, label1, node2.posNum, label2);
                        if (Math.abs(pairE) > maxInteraction) {
                            maxInteraction = Math.abs(pairE);
                        }
                    }
                }
                //Now we check if maxInteraction is greater than eCut for node1 and node2
                if (maxInteraction > eCut) {
                    interactionGraph[nodeNum1][nodeNum2] = true;
                    interactionGraph[nodeNum2][nodeNum1] = true;
                    countInteraction++;
                }
                possibleInteraction++;
            }
        }
        //System.out.println("Markov Random Field has " + countInteraction
//        +" pairs out of " + possibleInteraction + " possible pairs");
    }

    private ArrayList<MRFNode> getNeighbors(MRFNode node, boolean[][] interactionGraph) {
        ArrayList<MRFNode> neighbors = new ArrayList<>();
        for (MRFNode neighbor : this.nodeList) {
            //check if node is neighbor with node indexed by nodeIndex
            if (interactionGraph[node.posNum][neighbor.posNum]) {
                neighbors.add(neighbor);
            }
        }
        return neighbors;
    }

}
