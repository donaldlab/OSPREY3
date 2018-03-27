package edu.duke.cs.osprey.markstar.framework;

import java.util.PriorityQueue;

public class MARKStarNode implements Comparable<MARKStarNode> {
    private double upperBound;
    private double lowerBound;
    private PriorityQueue<MARKStarNode> children; // TODO: Pick appropriate data structure

    public MARKStarNode(){
        this.upperBound = Double.NaN;
        this.lowerBound = Double.NaN;
        this.children = new PriorityQueue<MARKStarNode>();
    }

    public MARKStarNode(double upperBound, double lowerBound, PriorityQueue<MARKStarNode> children){
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
        this.children = children;
    }

    public double getLowerBound(){
        return this.lowerBound;
    }

    public double getUpperBound(){
        return this.upperBound;
    }

    public PriorityQueue<MARKStarNode> getChildren(){
        return this.children;
    }

    public void updateBounds(double upperBound, double lowerBound){
        this.upperBound = upperBound;
        this.lowerBound = lowerBound;
    }

    @Override
    public int compareTo(MARKStarNode other){
        /**
         * Compares to another MARKStarNode.
         * The node with a larger difference between upper and lower bound is less.
         *
         * @param   other the {@code MARKStarNode} to compare
         * @return  1 if the other difference is smaller, -1 if the other difference
         *          is larger, and 0 otherwise.
         */
        if ((this.upperBound - this.lowerBound) < (other.upperBound - other.lowerBound)){
            return 1;
        }else if ((this.upperBound - this.lowerBound) > (other.upperBound - other.lowerBound)){
            return -1;
        }else {
            return 0;
        }
    }



}
