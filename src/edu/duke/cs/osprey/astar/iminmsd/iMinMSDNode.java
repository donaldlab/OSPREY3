/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.iminmsd;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTreeSuper;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * @author hmn5
 */
public class iMinMSDNode extends AStarNode {

    PruningMatrix[] pruneMat;
    ConfTreeSuper[] stateTrees;
    
    double boundGMEC;
    double unboundGMEC;
    
    public iMinMSDNode(int[] nodeAssignments, PruningMatrix[] pruneMat) {
        super(nodeAssignments, Double.NaN, false);
        this.pruneMat = pruneMat;
        boundGMEC = Double.NEGATIVE_INFINITY;
        unboundGMEC = Double.NEGATIVE_INFINITY;
    }
}
