package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.confspace.SearchProblem;

import java.util.Collection;
import java.util.PriorityQueue;
import java.util.Set;

public class MARKStarBound {

    /**
     * TODO: 1. Make MARKStarBounds use and update a queue.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */
    // We keep track of the root node for computing our K* bounds
    private MARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    private PriorityQueue<MARKStarNode> queue;
    private double epsilonBound = Double.POSITIVE_INFINITY;
    private boolean boundChanged = false;

    public MARKStarBound(SearchProblem problem){
        this.queue = new PriorityQueue<MARKStarNode>();
        rootNode = MARKStarNode.makeRoot(problem);
        queue.add(rootNode);
    }

    public double getBound(){
        if(boundChanged)
            recomputeBound(rootNode);
        return epsilonBound;
    }

    private void recomputeBound(MARKStarNode rootNode) {
    }

    public void tightenBound(){
        MARKStarNode nextNode = queue.poll();
        nextNode.expand();
        Collection<MARKStarNode> children = nextNode.getChildren();
        for(MARKStarNode child: children)
            child.computeBounds();
        boundChanged = true;

    }

}
