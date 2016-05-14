/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

/**
 *
 * @author mhall44
 */
public class AStarNode implements Comparable<AStarNode> {
    
    private int nodeAssignments[];//assignments (e.g. partial conformation) for node
    
    // partially-computed undefined energies (indexed by pos, rc)
    private double undefinedRCEnergies[][];
    
    private double score;//score (probably a lower bound on the energy)
    private double gscore;
    private double hscore;
    
    boolean scoreNeedsRefinement;

    
    //These are used in COMETS
    public double UB = Double.POSITIVE_INFINITY;//upper bound
    public int UBConf[] = null;//can have an upper bound on GMEC energy for this node's conf space
    //(and thus on the overall GMEC energy)
    //that is the energy of the conf denoted by UBConf
    
    
    
    //indicates the score needs to be refined (e.g. with EPIC continuous terms)
    //always false in simpler versions of A*
    public AStarNode(int[] nodeAssignments, boolean scoreNeedsRefinement) {
        this.nodeAssignments = nodeAssignments;
        this.undefinedRCEnergies = new double[this.nodeAssignments.length][];
        this.score = Double.NaN;
        this.gscore = Double.NaN;
        this.hscore = Double.NaN;
        this.scoreNeedsRefinement = scoreNeedsRefinement;
    }

    @Override
    public int compareTo(AStarNode other) {
        return Double.valueOf(score).compareTo(other.score);
    }

    public int[] getNodeAssignments() {
        return nodeAssignments;
    }
    
    public double[] getUndefinedRCEnergies(int pos) {
    	return undefinedRCEnergies[pos];
    }
    public void setUndefinedRCEnergies(int pos, double[] val) {
    	undefinedRCEnergies[pos] = val;
    }

    public void setScore(double score) {
        this.score = score;
    }
    public double getScore() {
        return score;
    }
    
    public double getGScore() {
    	return gscore;
    }
    public void setGScore(double val) {
    	gscore = val;
    }
    
    public double getHScore() {
    	return hscore;
    }
    public void setHScore(double val) {
    	hscore = val;
    }
    
    
    public boolean isFullyDefined(){
        //Assuming assignments greater than 0 denote fully defined positions,
        //determine if this node is fully defined or not
        for(int a : nodeAssignments){
            if(a<0)
                return false;
        }
        return true;
    }
    
    
}
