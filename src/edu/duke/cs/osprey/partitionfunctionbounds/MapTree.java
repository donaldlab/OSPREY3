/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.partitionfunctionbounds;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * @author hmn5
 */
public class MapTree extends ConfTree {

    int[] bestConf;
    double bestConfScore;

    public MapTree(SearchProblem sp) {
        super(sp);
    }

    public MapTree(EnergyMatrix emat, PruningMatrix pruneMat){
        super(emat, pruneMat);
    }

    @Override
    public boolean canPruneNode(AStarNode node) {
        if (node.getScore() > bestConfScore) {
            return true;
        }
        RCTuple tup = new RCTuple(node.getNodeAssignments());
        if (this.pruneMat.isPruned(tup)) {
            return true;
        }
        return false;
    }

}
