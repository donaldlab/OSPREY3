package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode;
import edu.duke.cs.osprey.sharkstar.SHARKStar;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.util.*;

public class SHARKStarNode implements ConfAStarNode {
    /**
     * Node in a Multi-sequence conformation tree
     *
     * TODO: Implement upper bound corrections
     * TODO: probably could optimize the memory usage
     */
    // Final variables
    public static final int Unassigned = -1;
    private static final boolean debug = MultiSequenceSHARKStarBound_refactor.debug;

    // Conformation variables
    private int[] assignments;

    // Score variables
    private double partialConfLB;                           // the pairwise-minimized free energy of the assigned residues
    private double partialConfUB;                           // the rigid pairwise free energy of the assigned residues
    private final Map<Sequence, Double> unassignedConfLB;   // the lower bound on the free energy of all full conformations
                                                                // compatible with this node (by sequence)
    private final Map<Sequence, Double> unassignedConfUB;   // the lower bound on the free energy of all full conformations
                                                                // compatible with this node (by sequence)
    private final Map<Sequence,Double> score;               // The ln of the partition function approximation error

    //TODO: Determine whether memory optimization requires deleting HOTCorrection, minE, and isMinimized
    private double HOTCorrectionLB;
    private double HOTCorrectionUB;
    private boolean isCorrectedLB;
    private boolean isCorrectedUB;
    private double minE;
    private boolean isMinimized;

    // Tree variables
    private final int level;
    private final SHARKStarNode parent;
    private final List<SHARKStarNode> children;

    // Debug variables
    private final List<String> history;


    public SHARKStarNode(int[] assignments, int level, SHARKStarNode parent){
        this.assignments = assignments;

        this.partialConfLB = Double.NaN;
        this.partialConfUB = Double.NaN;
        this.unassignedConfLB = new HashMap<>();
        this.unassignedConfUB = new HashMap<>();
        this.score = new HashMap<>();
        this.HOTCorrectionLB = 0.0;
        this.HOTCorrectionUB = 0.0;
        this.isCorrectedLB = false;
        this.isCorrectedUB = false;
        this.minE = Double.NaN;
        this.isMinimized = false;

        this.level = level;
        this.parent = parent;
        this.children = new ArrayList<>();

        if(debug){
            this.history = new ArrayList<>();
            if(this.parent != null) {
                history.add(String.format("Created node from parent %s", SimpleConfSpace.formatConfRCs(parent.assignments)));
            }else{
                history.add("Created root node.");
            }
        }
    }

    /**
     * Make a child node from this node
     */
    @Override
    public SHARKStarNode assign(int pos, int rc) {
        int[] newAssignments = new int[this.assignments.length];
        System.arraycopy(this.assignments, 0, newAssignments, 0, this.assignments.length);
        newAssignments[pos] = rc;
        SHARKStarNode child = new SHARKStarNode(newAssignments, this.level + 1, this);
        // Store the new node as a child
        this.children.add(child);
        return child;
    }


    @Override
    public void getConf(int[] conf) {
        // I don't know what this method is supposed to do
        throw new NotImplementedException();
    }

    public int[] getAssignments(){
        return this.assignments;
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

    @Override
    public double getGScore() {
        // Needed for calcOrder and calcDifferential
        return this.partialConfLB;
    }

    @Override
    public double getRigidGScore(){
        // needed for calcDifferential
        return this.partialConfUB;
    }

    @Override
    public void setGScore(double val) {
        System.err.println("ERROR: Don't try to set GScore directly. Instead, set the energy bounds.");
        throw new NotImplementedException();
    }

    @Override
    public double getHScore() {
        throw new NotImplementedException();
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
        return this.score.get(seq);
    }

    public void setScore(double score, Sequence seq){
        this.score.put(seq, score);
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

    public void setMinE(double val){
        this.minE = val;
    }

    public double getMinE(){
        return this.minE;
    }

    public void setIsMinimized(boolean val){
        this.isMinimized = val;
    }

    public boolean isMinimized(){
        return this.isMinimized;
    }

    public Map<Sequence, Double> getUnassignedConfLB(){
        return this.unassignedConfLB;
    }

    public Map<Sequence, Double> getUnassignedConfUB(){
        return this.unassignedConfUB;
    }

    public List<SHARKStarNode> getChildren(){
        return this.children;
    }

    public double getFreeEnergyLB(Sequence seq){
        if (this.isMinimized){
            return this.minE;
        }else if (this.isCorrectedLB){
            //TODO: This will change when I implement upperbound corrections
            return this.partialConfLB + this.HOTCorrectionLB + this.unassignedConfLB.get(seq);
            //return this.partialConfLB + this.unassignedConfLB.get(seq);
        }else{
            return this.partialConfLB + this.unassignedConfLB.get(seq);
        }
    }

    public double getFreeEnergyUB(Sequence seq){
        if (this.isMinimized){
            return this.minE;
        }else if (this.isCorrectedUB){
            //TODO: This will change when I implement upperbound corrections
            return this.partialConfUB + this.unassignedConfUB.get(seq) + this.HOTCorrectionUB;
        }else{
            return this.partialConfUB + this.unassignedConfUB.get(seq);
        }
    }

    public double getHOTCorrectionLB(){
        return this.HOTCorrectionLB;
    }

    public double getHOTCorrectionUB(){
        return this.HOTCorrectionLB;
    }

    public String confToString() {
        StringBuilder out = new StringBuilder("(");
        for (int assignment : assignments) {
            out.append(assignment).append(", ");
        }
        out.append(")");
        return out.toString();
    }

    public String toSeqString(Sequence seq){
        String out = String.format("%s -> [%.3f + %.3f, %.3f + %.3f]",
                confToString(), getPartialConfLB()+ getUnassignedConfLB(seq), HOTCorrectionLB, getPartialConfUB() + getUnassignedConfUB(seq), HOTCorrectionUB);
        if (isMinimized){
            out += String.format(" -> (minimized) %.3f", getMinE());
        }
        return out;
    }

    /**
     * Makes a SHARKStarNode generated during flexible precomputation compatible with a new conformation space
     *
     * @param permutation   The permutation array that maps the precomputed flexible positions to the positions in the new confSpace
     *                      so, permutation[i] gives the new ConfSpace index for residue i.
     * @param size  The size of the new confSpace
     */
    public void makeNodeCompatibleWithConfSpace(int[] permutation, int size){
        // first change the confSearch information

        // we copy over the new RCs based on permutation information
        int[] newAssignments = new int[size];
        Arrays.fill(newAssignments, -1);
        for (int i =0; i < this.assignments.length; i++){
            newAssignments[permutation[i]] = this.assignments[i];
        }
        // Now I'm going to be hacky and just copy over the assignments
        this.assignments = newAssignments;
    }

    public void setHOTCorrectionLB(double HOTCorrectionLB) {
        this.HOTCorrectionLB = HOTCorrectionLB;
        this.isCorrectedLB = true;
    }

    public void setHOTCorrectionUB(double HOTCorrectionUB) {
        this.HOTCorrectionUB = HOTCorrectionUB;
        this.isCorrectedUB = true;
    }

    public void unsetHOTCorrectionUB() {
        this.HOTCorrectionUB = 0;
        this.isCorrectedUB = false;
    }

    public Set<Sequence> getSequences(){
        Set<Sequence> intersect = this.unassignedConfLB.keySet();
        intersect.retainAll(this.unassignedConfUB.keySet());
        return intersect;
    }
}