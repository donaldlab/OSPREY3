/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.markovrandomfield;

import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.SuperRC;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.confspace.SearchProblemSuper;
import java.util.ArrayList;
/**
 *
 * @author hmn5
 */
public class MarkovRandomField {
    ArrayList<MRFNode> nodeList;
    boolean[][] interactionGraph;
    
    int numNodes;
    
    EnergyMatrix emat;
    PruningMatrix pruneMat;
    public MarkovRandomField(SearchProblemSuper searchProblem, double eCut){
        this.emat = searchProblem.emat;
        this.pruneMat = searchProblem.pruneMat;

        ConfSpaceSuper cSpace = searchProblem.confSpaceSuper;
        int numNodes = cSpace.numPos;
        
        //create nodeList
        for (int pos=0; pos<numNodes; pos++){
            MRFNode node = new MRFNode(pos, pruneMat.unprunedRCsAtPos(pos));
            nodeList.add(node);
        }
        this.numNodes = nodeList.size();
        
        //create interaction graph
        this.interactionGraph = createEnergyInteractionGraph(eCut);
        //create neighborList for each node
        for (MRFNode node : this.nodeList){
            node.neighborList = getNeighbors(node, this.interactionGraph);
        }
    }
    
    private boolean[][] createEnergyInteractionGraph(double eCut){
        boolean[][] interactionGraph = new boolean[numNodes][numNodes];
        int countInteraction = 0;
        int possibleInteraction = 0;
        //initialize all values to false
        for (int nodeNum1=0; nodeNum1<numNodes; nodeNum1++){
            for (int nodeNum2=0; nodeNum2<numNodes; nodeNum2++){
                interactionGraph[nodeNum1][nodeNum2] = false;
            }
        }
        //now get maxInteraction and check if it is greater than eCut
        for (int nodeNum1=0; nodeNum1<numNodes; nodeNum1++){
            for (int nodeNum2=0; nodeNum2<nodeNum1; nodeNum2++){
                double maxInteraction = 0.0;
                MRFNode node1 = this.nodeList.get(nodeNum1);
                MRFNode node2 = this.nodeList.get(nodeNum2);
                for (MRFLabel label1 : node1.labelList){
                    for (MRFLabel label2 : node2.labelList){
                        double pairE = emat.getPairwise(node1.nodeNum, label1.labelNum, node2.nodeNum, label2.labelNum);
                        if (Math.abs(pairE) > maxInteraction){
                            maxInteraction = Math.abs(pairE);
                        }
                    }
                }
                //Now we check if maxInteraction is greater than eCut for node1 and node2
                if (maxInteraction > eCut){
                    interactionGraph[nodeNum1][nodeNum2] = true;
                    interactionGraph[nodeNum2][nodeNum1] = true;
                    countInteraction++;
                }
                else{
                    possibleInteraction++;
                }
            }
        }
        System.out.println("Markov Random Field has "+countInteraction+
                            " pairs out of "+possibleInteraction+" possible pairs");
        return interactionGraph;
    }
    
    private ArrayList<MRFNode> getNeighbors(MRFNode node, boolean[][] interactionGraph){
        ArrayList<MRFNode> neighbors = new ArrayList<>();
        for (int nodeIndex=0;nodeIndex<this.numNodes;nodeIndex++){
            //check if node is neighbor with node indexed by nodeIndex
            if (interactionGraph[node.nodeNum][nodeIndex]){
                neighbors.add(this.nodeList.get(nodeIndex));
            }
        }
        return neighbors;
    }
    
}
