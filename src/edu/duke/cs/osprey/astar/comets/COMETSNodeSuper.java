/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * @author hmn5
 */
public class COMETSNodeSuper extends AStarNode {

    //score will be the lower bound on the objective function LME, for the node's sequence space
    PruningMatrix pruneMat[];//pruning matrix for each state

    //These things are only needed (and defined) for fully defined sequences
    ConfTreeSuper[] stateTrees = null;//trees searching conf space for each state
    double stateUB[] = null;//upper bounds on GMEC for each state

    public COMETSNodeSuper(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments, Double.NaN, false);//score not assigned yet, and doesn't need refinement
        this.pruneMat = pruneMat;
    }

    void expandConfTree() {
        //Pick one of the conformational search trees and expand it
        int firstSplittableLevel = Integer.MAX_VALUE;
        int stateToSplit = -1;

        for (int state = 0; state < stateTrees.length; state++) {
            ConfTreeSuper stateTree = stateTrees[state];

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
        ConfTreeSuper treeToSplit = stateTrees[stateToSplit];

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
