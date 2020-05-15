package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.sharkstar.LowerBoundException;
import edu.duke.cs.osprey.sharkstar.SHARKStar;
import edu.duke.cs.osprey.sharkstar.UpperBoundException;
import edu.duke.cs.osprey.tools.BigMath;

import java.math.BigDecimal;
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

    // User settings
    private final double boundTolerance = -1e-12;

    // variables
    private Map<Sequence, List<String>> minimizedNodes;

    public SHARKStarTreeDebugger(MathContext mathContext){
        this.mc = mathContext;
        this.bc = new BoltzmannCalculator(this.mc);

        this.minimizedNodes = new HashMap<>();
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
        // Initialize trackers
        if(this.minimizedNodes.containsKey(seq)){
            // reset minimized nodes
            this.minimizedNodes.remove(seq);
        }else{
            this.minimizedNodes.put(seq, new ArrayList<>());
        }
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
        enforceNoDuplicateChildren(children, seq);
        enforceBoundSanity(parent, seq);
        enforceUpperBoundCorrections(parent);
    }

    public void leafDebugChecks(SHARKStarNode leaf, Sequence seq){
        enforceBoundSanityLeaf(leaf, seq);
        if (leaf.isMinimized())
            minimizedNodes.get(seq).add(leaf.toSeqString(seq));
        enforceUpperBoundCorrections(leaf);
    }

    /**
     * Test to ensure that parent energy lower-bounds are less than child energy lower-bounds
     * @param parent    The parent SHARKStarNode
     * @param children  A list of child SHARKStarNodes
     * @param seq       The sequence for which to test energy bounds
     * @throws LowerBoundException  Exception indicating an issue with a lower bound
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
     * @throws UpperBoundException  Exception indicating an issue with an upper bound
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

    public void enforceNoDuplicateChildren(List<SHARKStarNode> children, Sequence seq){
        ArrayList<int[]> childAssignmentsList = children.stream()
                .map(SHARKStarNode::getAssignments)
                .collect(Collectors.toCollection(ArrayList::new));
        Set<int[]> childAssignmentsSet = new HashSet<>(childAssignmentsList);
        if (childAssignmentsList.size() != childAssignmentsSet.size()){
            children.stream().map(c -> c.toSeqString(seq)).forEach(System.out::println);
            throw new RuntimeException("Detected duplicate children!");
        }
    }

    public void enforceBoundSanity(SHARKStarNode node, Sequence seq){
        if (!(node.getFreeEnergyLB(seq) <= node.getFreeEnergyUB(seq))){
            System.err.println("ERROR: Insane bounds!");
            System.out.println(node.toSeqString(seq));
            throw new LowerBoundException();
        }
    }

    public void enforceBoundSanityLeaf(SHARKStarNode leaf, Sequence seq){
        if (leaf.isMinimized()) {
            if (!((leaf.getPartialConfLB() + leaf.getUnassignedConfLB(seq) + leaf.getHOTCorrectionLB() <= leaf.getMinE()) &&
                    (leaf.getMinE() <= leaf.getPartialConfUB() + leaf.getUnassignedConfUB(seq) + leaf.getHOTCorrectionUB()))) {
                System.err.println("ERROR: Insane bounds!");
                System.out.println(leaf.toSeqString(seq));
                throw new LowerBoundException();
            }
        }else{
            enforceBoundSanity(leaf,seq);
        }
    }

    public void setPrecomputedSequence(Sequence seq){
        this.precomputedSequence = seq;
    }

    public void printLeaves(Sequence seq){
        this.minimizedNodes.get(seq).forEach(System.out::println);
    }

    /**
     * Compute the Z bound by recursively moving through tree
     *
     * Probably only works for single-sequence trees
     * @param startNode     The node for which to compute Z bounds
     * @param seq           The sequence for which to compute Z bounds
     * @param energyMapper  A function that maps nodes to energies
     * @return              A partition function bound based on the energyMapper
     */
    public BigDecimal recursiveComputeZBound(SHARKStarNode startNode, Sequence seq, Function<SHARKStarNode, Double> energyMapper){
        List<SHARKStarNode> children = startNode.getChildren();
        if(children == null || children.size() == 0){
            if(startNode.getSequences().contains(seq)){
                return bc.calc(energyMapper.apply(startNode));
            }else{
                // This might be the wrong behavior for multi-sequence trees
                throw new RuntimeException(String.format("Reached dead node at %s",startNode.confToString()));
            }
        }else{
            BigMath sum = new BigMath(mc).set(0.0);
            for( SHARKStarNode child : children){
                sum.add(recursiveComputeZBound(child, seq, energyMapper));
            }
            return sum.get();
        }

    }

    public BigDecimal recursiveComputeZBound(Sequence seq, Function<SHARKStarNode, Double> energyMapper){
        return recursiveComputeZBound(this.rootNode, seq, energyMapper);
    }

    public void enforceUpperBoundCorrections(SHARKStarNode node){
        if(node.getHOTCorrectionUB()>0){
            throw new RuntimeException("Found an increasing upperbound correction");
        }
        if(node.getHOTCorrectionUB() != 0 && !MultiSequenceSHARKStarBound_refactor.doUpperBoundCorrections){
            throw new RuntimeException("Setting HOTCorrections when we shouldn't");
        }

    }
}
