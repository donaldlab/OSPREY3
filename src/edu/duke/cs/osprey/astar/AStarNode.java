/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

/**
 *
 * @author mhall44
 */
public class AStarNode implements Comparable {
    
    int nodeAssignments[];//assignments (e.g. partial conformation) for node
    
    double score;//score (probably a lower bound on the energy)
    
    boolean scoreNeedsRefinement;

    
    //indicates the score needs to be refined (e.g. with EPIC continuous terms)
    //always false in simpler versions of A*
    public AStarNode(int[] nodeAssignments, double score, boolean scoreNeedsRefinement) {
        this.nodeAssignments = nodeAssignments;
        this.score = score;
        this.scoreNeedsRefinement = scoreNeedsRefinement;
    }

    @Override
    public int compareTo(Object o) {
        AStarNode node2 = (AStarNode)o;//we can only compare to other AStarNodes, and expect no other cases
        return Double.valueOf(score).compareTo(node2.score);
    }
    
    
    
    
    
    
}
