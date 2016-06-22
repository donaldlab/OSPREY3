/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.seqkstar;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.partfunc.PartFuncTree;
import edu.duke.cs.osprey.partitionfunctionbounds.TRBP;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import java.util.ArrayList;

/**
 *
 * @author hmn5
 */
public class SequenceNode extends AStarNode {

    public PruningMatrix[] pruneMats;
    public PartFuncTree[] stateTrees;
    public double[] stateLBFreeEnergy;
    public double[] stateUBFreeEnergy;

    public SequenceNode(int[] aaNumPerPos, PruningMatrix pm) {
        super(aaNumPerPos, Double.NaN, false);
        this.pruneMats = new PruningMatrix[1];
        this.pruneMats[0] = pm;
    }

    public SequenceNode(int[] aaNumPerPos, PruningMatrix[] pruneMats) {
        super(aaNumPerPos, Double.NaN, false);
        this.pruneMats = pruneMats;
    }

    public void expandConfTree() {
        TRBP.setNumEdgeProbUpdates(3);
        for (int state = 0; state < stateTrees.length; state++) {
            AStarNode node = stateTrees[state].getQueue().poll();
            if (node != null) {
                ArrayList<AStarNode> children = stateTrees[state].getChildren(node);
                this.stateLBFreeEnergy[state] = -stateTrees[state].getCurrentUpperBoundLogZ();
                this.stateUBFreeEnergy[state] = -stateTrees[state].getCurrentLowerBoundLogZ();
                System.out.println("State: " + state + " " + this.stateTrees[state].effectiveEpsilon);
                for (AStarNode child : children) {
                    stateTrees[state].pq.add(child);
                }
            }
        }
    }
    
    public double getEffectiveEpsilon(int state){
        return stateTrees[state].effectiveEpsilon;
    }
}
