/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.partfunc;

import edu.duke.cs.osprey.astar.AStarNode;

/**
 *
 * @author hmn5
 */
public class partFuncNode extends AStarNode {
    
    //Lower bound on logZ under the node
    double lbLogZ;
    //Upper bound on logZ under the node
    double ubLogZ;
    
    public double[][] edgeProbabilities;
    
    public partFuncNode(int[] anodeAssignments, double score, boolean scoreNeedsRefinement){
        super(anodeAssignments, score, scoreNeedsRefinement);
    }
    public partFuncNode(int[] anodeAssignments){
        super (anodeAssignments, Double.NaN, false);
    }
    
    public void setLowerBoundLogZ(double albLogZ){
        this.lbLogZ = albLogZ;
    }
    
    public void setUpperBoundLogZ(double aubLogZ){
        this.ubLogZ = aubLogZ;
    }
    
    public double getLowerBoundLogZ(){
        return this.lbLogZ;
    }
    
    public double getUpperBoundLogZ(){
        return this.ubLogZ;
    }
}
