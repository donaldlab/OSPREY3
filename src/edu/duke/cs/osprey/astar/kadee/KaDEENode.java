/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.kadee;

import edu.duke.cs.osprey.astar.AStarNode;
import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.pruning.PruningMatrix;
/**
 *The class implements a node for the KaDEETree
 * 
 * @author hmn5
 */
public class KaDEENode extends AStarNode{
    // The node score will be a lower bound on the change in free energy
    
    PruningMatrix[] pruneMat;//pruning matrxi for each state
    
    //When a sequence node is fully defined, it has confTree
    ConfTree[] stateTrees;
    
    public KaDEENode(int[] nodeAssignments, PruningMatrix[] pruneMat){
        super(nodeAssignments, Double.NaN, false);
        this.pruneMat = pruneMat;
    }
    
}
