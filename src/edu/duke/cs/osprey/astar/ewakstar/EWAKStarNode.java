package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.FullAStarNode;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.ArrayList;
import java.util.PriorityQueue;

public class EWAKStarNode extends FullAStarNode{
    //score will be the lower bound on the objective function LME, for the node's sequence space

    PruningMatrix pruneMat[];//pruning matrix for each state

    //These things are only needed (and defined) for fully defined sequences
    ConfTree<FullAStarNode>[] stateTrees = null;//trees searching conf space for each state
    double stateUB[] = null;//upper bounds on GMEC for each state



    public EWAKStarNode(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments);//score not assigned yet, and doesn't need refinement
        this.pruneMat = pruneMat;
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
