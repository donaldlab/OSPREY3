/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class SequenceNode extends AStarNode {

    public PruningMatrix pruneMat;
    public PartFuncTree tree;
    public double lowerBoundFreeEnergy;
    public double upperBoundFreeEnergy;

    public SequenceNode(int[] aaNumPerPos, PruningMatrix pm) {
        super(aaNumPerPos, Double.NaN, false);
        this.pruneMat = pm;
    }

    public void expandConfTree() {
        ArrayList<AStarNode> children = tree.getChildren(tree.getQueue().poll());
        this.lowerBoundFreeEnergy = -tree.getCurrentUpperBoundLogZ();
        this.upperBoundFreeEnergy = -tree.getCurrentLowerBoundLogZ();
        for (AStarNode child : children) {
            tree.pq.add(child);
        }
    }
}
