package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.*;

public class SHARKStarNode implements ConfAStarNode {
    /**
     * Node in a Multi-sequence conformation tree
     *
     * SHARK* explores the conformation tree in order of free energy approximation error (equivalently, partition
     * function approximation error). This means that the 'GScore' for a SHARK*Node is actually the difference between
     * the upper and lower bounds on the partial conformation energy. Similarly, the 'HScore' is the difference between
     * the upper and lower bounds on the unassigned conformation free energy. In essence, the scores for these nodes are
     * measures of *error*, not energy.
     *
     * TODO: Implement upper bound corrections
     * TODO: probably could optimize the memory usage
     */
    // Final variables
    private static final int Unassigned = -1;
    private static final boolean debug = true;

    // Conformation variables
    private int[] assignments;

    // Score variables
    private double partialConfLB; // the pairwise-minimized free energy of the assigned residues
    private double partialConfUB; // the rigid pairwise free energy of the assigned residues
    private Map<Sequence, Double> unassignedConfLB; // the lower bound on the free energy of all full conformations
                                                    // compatible with this node (by sequence)
    private Map<Sequence, Double> unassignedConfUB; // the lower bound on the free energy of all full conformations
                                                    // compatible with this node (by sequence)
    //TODO: Determine whether memory optimization requires deleting HOTCorrection, minE, and isMinimized
    private double HOTCorrection;
    private boolean isCorrected;
    private double minE;
    private boolean isMinimized;

    // Tree variables
    private int level;
    private SHARKStarNode parent;
    private Map<Sequence, List<SHARKStarNode>> children;

    // Debug variables
    private List<String> history;


    public SHARKStarNode(int[] assignments, int level, SHARKStarNode parent){
        this.assignments = assignments;

        this.partialConfLB = Double.NaN;
        this.partialConfUB = Double.NaN;
        this.unassignedConfLB = new HashMap<>();
        this.unassignedConfUB = new HashMap<>();
        this.HOTCorrection = 0.0;
        this.isCorrected = false;
        this.minE = Double.NaN;
        this.isMinimized = false;

        this.level = level;
        this.parent = parent;
        this.children = new HashMap<>();

        if(debug){
            this.history = new ArrayList<>();
            history.add(String.format("Created node from parent %s", SimpleConfSpace.formatConfRCs(parent.assignments)));
        }
    }

    /**
     * Make a child node from this node
     */
    @Override
    public SHARKStarNode assign(int pos, int rc) {
        System.err.println("WARNING: Making child without sequence. This will break the tree links.");

        int[] newAssignments = new int[this.assignments.length];
        System.arraycopy(this.assignments, 0, newAssignments, 0, this.assignments.length);
        newAssignments[pos] = rc;
        return new SHARKStarNode(newAssignments, this.level + 1, this);
    }

    /**
     * Make a child node from this node and assign it to this node's children list
     * TODO: Determine whether this method is useful
     */
    public SHARKStarNode assign(int pos, int rc, Sequence seq) {
        // Make the new node
        int[] newAssignments = new int[this.assignments.length];
        System.arraycopy(this.assignments, 0, newAssignments, 0, this.assignments.length);
        newAssignments[pos] = rc;
        SHARKStarNode child = new SHARKStarNode(newAssignments, this.level + 1, this);
        // Store the new node as a child
        this.children.get(seq).add(child);
        return child;
    }

    @Override
    public void getConf(int[] conf) {
        // I don't know what this method is supposed to do
        throw new NotImplementedException();
    }

    /**
     * Tell the confIndex about the information in this node
     */
    @Override
    public void index(ConfIndex index) {
        index.numDefined = 0;
        index.numUndefined = 0;
        for (int pos = 0; pos < assignments.length; pos++) {
            int rc = assignments[pos];
            if (rc == Unassigned) {
                index.undefinedPos[index.numUndefined] = pos;
                index.numUndefined++;
            } else {
                index.definedPos[index.numDefined] = pos;
                index.definedRCs[index.numDefined] = assignments[pos];
                index.numDefined++;
            }
        }
        index.node = this;
    }

    /**
     * Returns the free energy approximation error in the partial conformation
     */
    @Override
    public double getGScore() {
        if (this.isMinimized){
            return 0.0;
        }else if (this.isCorrected){
            //TODO: This will change when I implement upperbound corrections
            return this.partialConfUB - this.partialConfLB + this.HOTCorrection;
        }else{
            return this.partialConfUB - this.partialConfLB;
        }
    }

    @Override
    public void setGScore(double val) {
        System.err.println("ERROR: Don't try to set GScore directly. Instead, set the energy bounds.");
        throw new NotImplementedException();
    }

    @Override
    public double getHScore() {
        System.out.println("WARNING: Getting HScore of a Multi-sequence node without specifying sequence. This probably isn't what you want to do.");
        throw new NotImplementedException();
    }

    /**
     * Returns the free energy approximation error in the unassigned conformation
     */
    public double getHScore(Sequence seq) {
        if (this.isMinimized){
            return 0.0;
        }else if (this.isCorrected){
            //TODO: This will change when I implement upperbound corrections
            return this.unassignedConfUB.get(seq) - this.unassignedConfLB.get(seq);
        }else{
            return this.unassignedConfUB.get(seq) - this.unassignedConfLB.get(seq);
        }
    }

    @Override
    public void setHScore(double val) {
        System.err.println("ERROR: Don't try to set HScore directly. Instead, set the energy bounds.");
        throw new NotImplementedException();
    }

    @Override
    public double getScore(){
        System.out.println("WARNING: Getting Score of a Multi-sequence node without specifying sequence. This probably isn't what you want to do.");
        throw new NotImplementedException();
    }

    public double getScore(Sequence seq){
        return getGScore() + getHScore(seq);
    }

    @Override
    public int getLevel() {
        return this.level;
    }

    public void setPartialConfLB(double val){
        if(debug)
            this.history.add(String.format("Setting partialConfLB <- %.3f", val));
        this.partialConfLB = val;
    }

    public double getPartialConfLB(){
        return this.partialConfLB;
    }

    public void setPartialConfUB(double val){
        if(debug)
            this.history.add(String.format("Setting partialConfUB <- %.3f", val));
        this.partialConfUB = val;
    }

    public double getPartialConfUB(){
        return this.partialConfUB;
    }

    public void setUnassignedConfLB(double val, Sequence seq){
        if(debug)
            this.history.add(String.format("Setting unassignedConfLB( %s ) <- %.3f", seq, val));
            this.unassignedConfLB.put(seq, val);
    }

    public double getUnassignedConfLB(Sequence seq){
        return this.unassignedConfLB.get(seq);
    }

    public void setUnassignedConfUB(double val, Sequence seq){
        if(debug)
            this.history.add(String.format("Setting unassignedConfUB( %s ) <- %.3f", seq, val));
        this.unassignedConfUB.put(seq, val);
    }

    public double getUnassignedConfUB(Sequence seq){
        return this.unassignedConfUB.get(seq);
    }
}
