/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 * The class implements a node for the KaDEETree
 *
 * @author hmn5
 */
public class KaDEENode extends AStarNode {
    // The node score will be a lower bound on the change in free energy

    PruningMatrix[] pruneMat;//pruning matrxi for each state

    //When a sequence node is fully defined, it has confTree
    ConfTree[] stateTrees;
    double[] stateUB = null;//upper bounds on GMEC for each state

    public boolean scoreSet = false;
    
    public KaDEENode(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments, Double.NaN, false);
        this.pruneMat = pruneMat;
    }

    void expandConfTree() {
        //Pick one of the conformational search trees and expand it
        int firstSplittableLevel = Integer.MAX_VALUE;
        int stateToSplit = -1;

        for (int state = 0; state < stateTrees.length; state++) {
            ConfTree stateTree = stateTrees[state];

            if (stateTree != null) {
                AStarNode curBestNode = stateTree.getQueue().peek();
                if (!curBestNode.isFullyDefined()) {
                    int stateFirstSplittableLevel = stateTree.nextLevelToExpand(curBestNode.getNodeAssignments());
                    if (stateFirstSplittableLevel < firstSplittableLevel) {
                        //tree has a splittable level, and it's the lowest so far
                        stateToSplit = state;
                        firstSplittableLevel = stateFirstSplittableLevel;
                    }
                }
            }
        }

        //if(stateToSplit==-1)//all trees fully defined!
        //    return;
        ConfTree treeToSplit = stateTrees[stateToSplit];

        PriorityQueue<AStarNode> expansion = treeToSplit.getQueue();

        AStarNode bestNode = expansion.poll();
        ArrayList<AStarNode> children = treeToSplit.getChildren(bestNode);

        for (AStarNode child : children) {
            child.UB = Double.NaN;//update UB later if needed
            child.UBConf = bestNode.UBConf;//store this as an initial value
            expansion.add(child);
        }
    }
    

}
