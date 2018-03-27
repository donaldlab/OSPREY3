package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.markstar.MARKStar;

import java.util.PriorityQueue;

public class MARKStarBound {
    private MARKStarNode curMNode;
    private PriorityQueue<MARKStarNode> queue;

    public MARKStarBound(){
        this.curMNode = new MARKStarNode();
        this.queue = new PriorityQueue<MARKStarNode>();
    }

    public MARKStarBound(MARKStarNode root){
        this.curMNode = root;
        this.queue = new PriorityQueue<MARKStarNode>();
    }

    public double getBound(){
        //TODO: implement method
        /*
        Pseudocode:
        computeUpperBound()
        computeLowerBound()
        return uBound -lBound
         */
        return Double.NaN;
    }

    public void tightenBound(){
        //TODO: implement method
        /*
        Pseudocode:
        curMNode = queue.pop();
        expand(curMNode)
        curMNode.updateBounds();
         */
    }
}
