package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.sharkstar.LowerBoundException;
import edu.duke.cs.osprey.sharkstar.SHARKStar;
import edu.duke.cs.osprey.sharkstar.UpperBoundException;

import java.math.MathContext;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class SHARKStarTreeDebugger {
    /**
     * Class to check correctness of a tree or subtree of SHARKStarNode instances
     *
     */

    private final MathContext mc;
    private final BoltzmannCalculator bc;
    private SHARKStarNode rootNode;
    private Sequence precomputedSequence;

    private final double boundTolerance = -1e-12;

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
        // skip the root node because we have issues with weird sequences and such
        for( SHARKStarNode child : rootNode.getChildren()){
            travelNode(child, seq);
        }
    }

    /**
     * Run debug checks, then recurse to children
     *
     * @param node The node to check
     * @param seq The sequence to look at
     */
    public void travelNode(SHARKStarNode node, Sequence seq){
        boolean nodeContainsSeq = node.getSequences().contains(seq);
        // Get the children
        List<SHARKStarNode> children = node.getChildren();
        // if there are multiple children
        if(children.size() > 0) {
            // run debug checks on the children if the node has an entry for seq
            if (nodeContainsSeq)
                debugChecks(node, children, seq);
            // recurse
            for (SHARKStarNode child : children) {
                travelNode(child, seq);
            }
        }
        // If we are at a leaf, run different tests
        if (nodeContainsSeq)
            leafDebugChecks(node, seq);
    }

    public void debugChecks(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq){
        enforceDecreasingUpperBounds(parent, children, seq);
        enforceIncreasingLowerBounds(parent, children, seq);
        enforceNoDuplicateChildren(children);
    }

    public void leafDebugChecks(SHARKStarNode leaf, Sequence seq){
    }

    /**
     * Test to ensure that parent energy lower-bounds are less than child energy lower-bounds
     * @param parent    The parent SHARKStarNode
     * @param children  A list of child SHARKStarNodes
     * @param seq       The sequence for which to test energy bounds
     * @throws LowerBoundException
     */
    public void enforceIncreasingLowerBounds(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq) throws LowerBoundException{
        // compute the child energies for sequence children
        double childE;
        if(children.size() > 1){
            ArrayList<Double> childEnergies = children.stream()
                    .filter(n -> n.getSequences().contains(seq))
                    .map(n -> n.getFreeEnergyLB(seq))
                    .collect(Collectors.toCollection(ArrayList::new));
            childE= bc.logSumExp(childEnergies);
        }else{
            childE= children.get(0).getFreeEnergyLB(seq);
        }

        // Compute the parent energy
        double parentE = parent.getFreeEnergyLB(seq);

        // check to make sure the bounds are decreasing
        if (!(childE - parentE > boundTolerance)){ // parentE < childE
            //Print info
            System.err.println(String.format("ERROR: Lower bounds are decreasing from parent to children! %.16f > %.16f", parentE, childE));
            System.out.println(String.format("Difference of %2.3e",childE - parentE));
            List<SHARKStarNode> parentList = new ArrayList<>();
            parentList.add(parent);
            Stream.of(parentList, children)
                    .flatMap(Collection::stream)
                    .filter(n -> n.getSequences().contains(seq))
                    .forEach(e -> System.out.println(e.toSeqString(seq)));
            throw new LowerBoundException();
        }
    }

    /**
     * Test to ensure that parent energy upper-bounds are greater than child energy upper-bounds
     * @param parent    The parent SHARKStarNode
     * @param children  A list of child SHARKStarNodes
     * @param seq       The sequence for which to test energy bounds
     * @throws UpperBoundException
     */
    public void enforceDecreasingUpperBounds(SHARKStarNode parent, List<SHARKStarNode> children, Sequence seq) throws UpperBoundException{
        // compute the child energies for sequence children
        double childE;
        if(children.size() > 1){
            ArrayList<Double> childEnergies = children.stream()
                    .filter(n -> n.getSequences().contains(seq))
                    .map(n -> n.getFreeEnergyUB(seq))
                    .collect(Collectors.toCollection(ArrayList::new));
            childE= bc.logSumExp(childEnergies);
        }else{
            childE= children.get(0).getFreeEnergyUB(seq);
        }

        // Compute the parent energy
        double parentE = parent.getFreeEnergyUB(seq);

        // check to make sure the bounds are decreasing
        if (!(parentE - childE > boundTolerance)){ // parentE > childE
            //Print info
            System.err.println(String.format("ERROR: Upper bounds are increasing from parent to children! %.16f < %.16f", parentE, childE));
            System.out.println(String.format("Difference of %2.3e",parentE - childE));
            List<SHARKStarNode> parentList = new ArrayList<>();
            parentList.add(parent);
            Stream.of(parentList, children)
                    .flatMap(Collection::stream)
                    .filter(n -> n.getSequences().contains(seq))
                    .forEach(e -> System.out.println(e.toSeqString(seq)));
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

    public void setPrecomputedSequence(Sequence seq){
        this.precomputedSequence = seq;
    }
}
