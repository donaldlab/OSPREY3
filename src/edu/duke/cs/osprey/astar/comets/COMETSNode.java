/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import java.util.ArrayList;
import java.util.PriorityQueue;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * Node for sequence-based A* (COMETS) tree
 * Has to deal with conformational search, so significantly more complex than standard A* node
 * 
 * @author mhall44
 */
public class COMETSNode extends FullAStarNode {
    
    
    //score will be the lower bound on the objective function LME, for the node's sequence space
    
    PruningMatrix pruneMat[];//pruning matrix for each state
        
    //These things are only needed (and defined) for fully defined sequences
    ConfTree<FullAStarNode>[] stateTrees = null;//trees searching conf space for each state
    double stateUB[] = null;//upper bounds on GMEC for each state

    
    
    public COMETSNode(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments);//score not assigned yet, and doesn't need refinement
        this.pruneMat = pruneMat;
    }
    
    public COMETSNode(COMETSNode cn){
        //copy constructor.  COMETS-specific things are big so let's not copy them unnecessarily
        super(cn);
        pruneMat = cn.pruneMat;
        stateTrees = cn.stateTrees;
        stateUB = cn.stateUB;
    }
    
    void expandConfTree(){
        //Pick one of the conformational search trees and expand it
        int firstSplittableLevel = Integer.MAX_VALUE;
        int stateToSplit = -1;

        
        for(int state=0; state<stateTrees.length; state++){
            ConfTree<FullAStarNode> stateTree = stateTrees[state];

            if(stateTree!=null){
                FullAStarNode curBestNode = stateTree.getQueue().peek();
                if( ! curBestNode.isFullyDefined() ){
                    int stateFirstSplittableLevel = stateTree.nextLevelToExpand( curBestNode );
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

        ConfTree<FullAStarNode> treeToSplit = stateTrees[stateToSplit];
        
        PriorityQueue<FullAStarNode> expansion = treeToSplit.getQueue();
        
        FullAStarNode bestNode = (FullAStarNode)expansion.poll();
        ArrayList<FullAStarNode> children = treeToSplit.getChildren(bestNode);
        
        for(FullAStarNode child : children){
            child.UB = Double.NaN;//update UB later if needed
            child.UBConf = bestNode.UBConf;//store this as an initial value
            expansion.add(child);
        }
    }
    
    
    
    
    
}
