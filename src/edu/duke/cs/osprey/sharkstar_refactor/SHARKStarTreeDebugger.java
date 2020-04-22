package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.sharkstar.LowerBoundException;
import edu.duke.cs.osprey.sharkstar.SHARKStar;
import edu.duke.cs.osprey.sharkstar.UpperBoundException;

import java.math.MathContext;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SHARKStarTreeDebugger {
    /**
     * Class to check correctness of a tree or subtree of SHARKStarNode instances
     *
     */

    private final MathContext mc;
    private final BoltzmannCalculator bc;
    private SHARKStarNode rootNode;

    public SHARKStarTreeDebugger(MathContext mathContext){
        this.mc = mathContext;
        this.bc = new BoltzmannCalculator(this.mc);
    }

    public void setRootNode(SHARKStarNode node){
        this.rootNode = node;
    }

    /**
     * Recurse down from the root node, running debug checks along the way
     *
     * @param seq The sequence to look at
     */
    public void travelTree(Sequence seq){
        travelNode(this.rootNode, seq);
    }

    /**
     * Run debug checks, then recurse to children
     *
     * @param node The node to check
     * @param seq The sequence to look at
     */
    public void travelNode(SHARKStarNode node, Sequence seq){
        // Get the children
        //TODO: Make this get children compatible with sequence
        List<SHARKStarNode> children = node.getChildren();
        // if there are children
        if(children != null) {
            // run debug checks on the children
            debugChecks(node, children, seq);
            // recurse
            for (SHARKStarNode child : children) {
                travelNode(child, seq);
            }
        }
        // If we are at a leaf, run different tests
        leafDebugChecks(node, seq);
    }

    public void debugChecks(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq){
        enforceConsistentUpperBounds(parent, children, seq);
        enforceConsistentLowerBounds(parent, children, seq);
        enforceNoDuplicateChildren(children);
    }

    public void leafDebugChecks(SHARKStarNode leaf, Sequence seq){
    }

    public void enforceConsistentUpperBounds(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq){
        ArrayList<Double> childEnergies = children.stream()
                .map(n -> n.getFreeEnergyLB(seq))
                .collect(Collectors.toCollection(ArrayList::new));
        if (parent.getFreeEnergyLB(seq) > bc.logSumExp(childEnergies)){
            throw new LowerBoundException();
        }
    }

    public void enforceConsistentLowerBounds(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq){
        ArrayList<Double> childEnergies = children.stream()
                .map(n -> n.getFreeEnergyUB(seq))
                .collect(Collectors.toCollection(ArrayList::new));
        if (parent.getFreeEnergyUB(seq) < bc.logSumExp(childEnergies)){
            throw new UpperBoundException();
        }
    }

    public void enforceNoDuplicateChildren(List<SHARKStarNode> children){
        ArrayList<int[]> childAssignmentsList = children.stream()
                .map(SHARKStarNode::getAssignments)
                .collect(Collectors.toCollection(ArrayList::new));
        Set<int[]> childAssignmentsSet = new HashSet<>(childAssignmentsList);
        if (childAssignmentsList.size() != childAssignmentsSet.size()){
            throw new RuntimeException("Detected duplicate children!");
        }
    }
}
