/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 *
 * @author hmn5
 */
public class SequenceNode extends AStarNode {
    public PruningMatrix pruneMat;
    
    public SequenceNode(int[] aaNumPerPos, PruningMatrix pm){
        super(aaNumPerPos, Double.NaN, false);
        this.pruneMat = pm;
    }
}
