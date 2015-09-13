/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;
import java.util.PriorityQueue;

/**
 *
 * Node for sequence-based A* (COMETS) tree
 * Has to deal with conformational search, so significantly more complex than standard A* node
 * 
 * @author mhall44
 */
public class COMETSNode extends AStarNode {
    
    
    //score will be the lower bound on the objective function LME, for the node's sequence space
    
    PruningMatrix pruneMat[];//pruning matrix for each state
        
    //These things are only needed (and defined) for fully defined sequences
    ConfTree[] stateTrees = null;//trees searching conf space for each state
    double stateUB[] = null;//upper bounds on GMEC for each state

    
    
    public COMETSNode(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments, Double.NaN, false);//score not assigned yet, and doesn't need refinement
        this.pruneMat = pruneMat;
    }
    
    
    
    void expandConfTree(){
        //Pick one of the conformational search trees and expand it
        int firstSplittableLevel = Integer.MAX_VALUE;
        int stateToSplit = -1;

        
        for(int state=0; state<stateTrees.length; state++){
            ConfTree stateTree = stateTrees[state];

            if(stateTree!=null){
                AStarNode curBestNode = stateTree.getQueue().peek();
                if( ! curBestNode.isFullyDefined() ){
                    int stateFirstSplittableLevel = stateTree.nextLevelToExpand( curBestNode.getNodeAssignments() );
                    if( stateFirstSplittableLevel < firstSplittableLevel ){
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
        
        for(AStarNode child : children){
            child.UB = Double.NaN;//update UB later if needed
            child.UBConf = bestNode.UBConf;//store this as an initial value
            expansion.add(child);
        }
    }
    
    
    
    
    
}
