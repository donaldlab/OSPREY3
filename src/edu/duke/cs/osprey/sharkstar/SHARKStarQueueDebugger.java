package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.ObjectPool;
import jdk.internal.net.http.common.Pair;

import java.math.BigDecimal;
import java.util.*;
import java.util.stream.Collectors;

public class SHARKStarQueueDebugger {
    /**
     * A class to debug a SHARKStarQueue_refactor object
     */

    private static BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    public SHARKStarQueueDebugger(){

    }

    /**
     * Checks for duplicates in the SHARKStarQueue
     *
     * @param queue The queue to check
     * @return true if there are duplicates, false otherwise
     */
    public static boolean hasDuplicateNodes(Collection<MultiSequenceSHARKStarNode> queue){
        return queue.parallelStream().map( node -> new AbstractMap.SimpleEntry<>(node, 1)) // map nodes to nodes with freq counter
                .collect(Collectors.groupingBy(e -> e.getKey().assignments, Collectors.counting())) //collect into map with freqs
                .values().stream().map(e -> e>1).reduce(false, (a,b) -> a || b); // check whether any frequencies are greater than 1
    }

    /**
     * Prints the number of nodes at each level in queue
     *
     * @param queue The queue to check
     */
    public static void printLevelBreakdown(Collection<MultiSequenceSHARKStarNode> queue){
        queue.parallelStream().map(node -> new AbstractMap.SimpleEntry<>(node.getLevel(), 1)) // map nodes to level and freq counter
                .collect(Collectors.groupingBy(AbstractMap.SimpleEntry::getKey, Collectors.counting()))
                .forEach((k, v) -> System.out.println(String.format("Level %d: %d nodes", k, v)));
    }

    /**
     * Prints the node information for sequence for each node at level
     *
     * @param queue The queue to check
     * @param level The level at which to print nodes
     * @param seq   The sequence for which to print information
     */
    public static void printNodesAtLevel(Collection<MultiSequenceSHARKStarNode> queue, int level, Sequence seq){
        System.out.println(String.format("Printing all nodes at level %d:", level));
        List<MultiSequenceSHARKStarNode> nodes = getNodesAtLevel(queue, level);
        nodes.sort(Comparator.comparing(a -> a.getErrorBound(seq)));
        nodes.forEach(node -> System.out.println(String.format("%s: [%.3f, %.3f], [%1.3e, %1.3e]",
                node.confToString(),
                node.getConfLowerBound(seq),
                node.getConfUpperBound(seq),
                bc.calc(node.getConfUpperBound(seq)),
                bc.calc(node.getConfLowerBound(seq))
        )));
    }

    /**
     * Returns the nodes at level
     *
     * @param queue The queue to check
     * @param level The level at which to print nodes
     */
    public static List<MultiSequenceSHARKStarNode> getNodesAtLevel(Collection<MultiSequenceSHARKStarNode> queue, int level) {
        return queue.parallelStream().filter(node -> node.getLevel() == level)
                .collect(Collectors.toList());
    }

    public static void compareNodes(List<MultiSequenceSHARKStarNode> list1, List<MultiSequenceSHARKStarNode> list2, Sequence seq){
        if(hasDuplicateNodes(list1))
            throw new RuntimeException("Error: first argument has duplicates");
        else if (hasDuplicateNodes(list2))
            throw new RuntimeException("Error: second argument has duplicates");
        Map<int[], MultiSequenceSHARKStarNode> firstAssignments = new HashMap<>();
        Map<int[], List<MultiSequenceSHARKStarNode>> secondAssignments = new HashMap<>();
        list1.forEach((a) -> firstAssignments.put(a.assignments, a));

        for(MultiSequenceSHARKStarNode node : list2){
            boolean foundParent=false;
            for(int [] ass : firstAssignments.keySet()){
                if(isPotentialChildOf(ass, node.assignments)){
                    secondAssignments.putIfAbsent(ass, new ArrayList<>());
                    secondAssignments.get(ass).add(node);
                    foundParent=true;
                }
            }
            if(!foundParent){
                System.err.println(String.format("Did not find parent for %s", node.confToString()));
            }
        }
        //list2.forEach((a) -> secondAssignments.put(a.assignments, a));

        list1.sort(Comparator.comparing((a) -> a.getErrorBound(seq)));
        for (MultiSequenceSHARKStarNode node : list1){
            if(secondAssignments.containsKey(node.assignments)) {
                List<MultiSequenceSHARKStarNode> children = secondAssignments.get(node.assignments);
                children.sort(Comparator.comparing(MultiSequenceSHARKStarNode::getLevel));
                BigDecimal childError = children.stream().map((n) -> bc.calc(n.getErrorBound(seq))).
                        reduce(BigDecimal.ZERO, BigDecimal::add);
                BigDecimal childLB = children.stream().map((n) -> bc.calc(n.getConfUpperBound(seq))).
                        reduce(BigDecimal.ZERO, BigDecimal::add);
                BigDecimal childUB = children.stream().map((n) -> bc.calc(n.getConfLowerBound(seq))).
                        reduce(BigDecimal.ZERO, BigDecimal::add);
                BigDecimal deltaError = bc.calc(node.getErrorBound(seq)).subtract(childError);
                if(deltaError.compareTo(BigDecimal.ZERO) < 0){
                    System.err.println("Bounds increasing!?");
                    System.err.println(String.format("Delta error %1.3e: [%1.3e, %1.3e] -> [%1.3e, %1.3e]",
                            deltaError,
                            bc.calc(node.getConfUpperBound(seq)),
                            bc.calc(node.getConfLowerBound(seq)),
                            childLB,
                            childUB
                    ));
                }
                if(deltaError.compareTo(BigDecimal.ONE) > 0) {
                    System.out.println(String.format("Delta error %1.3e: [%1.3e, %1.3e] -> [%1.3e, %1.3e]",
                            deltaError,
                            bc.calc(node.getConfUpperBound(seq)),
                            bc.calc(node.getConfLowerBound(seq)),
                            childLB,
                            childUB
                    ));
                    System.out.println(String.format("%s: [%.3f, %.3f] --> [%1.3e, %1.3e]",
                            node.confToString(),
                            node.getConfLowerBound(seq),
                            node.getConfUpperBound(seq),
                            bc.calc(node.getConfUpperBound(seq)),
                            bc.calc(node.getConfLowerBound(seq))
                            ));
                    for (MultiSequenceSHARKStarNode child : children){
                        System.out.println(String.format("\t%s: [%.3f, %.3f] --> [%1.3e, %1.3e]",
                                child.confToString(),
                                child.getConfLowerBound(seq),
                                child.getConfUpperBound(seq),
                                bc.calc(child.getConfUpperBound(seq)),
                                bc.calc(child.getConfLowerBound(seq))
                        ));
                    }
                }
            }
        }
    }

    public static boolean isValidFringeForSeq(Collection<MultiSequenceSHARKStarNode> fringe, Sequence seq, RCs seqRCs){
        Map<Integer, List<MultiSequenceSHARKStarNode>> nodesByLevel = fringe.stream().collect(Collectors.groupingBy(MultiSequenceSHARKStarNode::getLevel));
        // iterate over the positions
        for (int level = seqRCs.getNumPos(); level > 0; level--){
            // skip positions for which there are no nodes
            if(!nodesByLevel.containsKey(level))
                continue;
            int posIndex = nodesByLevel.get(level).get(0).pos; // get the pos by looking at the first node
            // group the children by their parent assignments
            Map<List<Integer>, List<MultiSequenceSHARKStarNode>> nodesByParentAssignments = nodesByLevel.get(level).stream()
                    .collect(Collectors.groupingBy((n) ->{
                        List<Integer> parentAssignments = Arrays.stream(n.assignments).boxed().collect(Collectors.toList());
                        parentAssignments.set(posIndex, MultiSequenceSHARKStarNode.Unassigned);
                        return parentAssignments;
                    }));

            // check to make sure that we have the full set of children for each parent
            for(List<Integer> parentAssignments : nodesByParentAssignments.keySet()){
                List<MultiSequenceSHARKStarNode> children = nodesByParentAssignments.get(parentAssignments);
                if(children.size() < seqRCs.getNum(posIndex))
                    return false;
                // Create the parent from the children and add it to the nodesByLevel
                MultiSequenceSHARKStarNode parent = new MultiSequenceSHARKStarNode(seqRCs.getNumPos(), level - 1);
                int[] arrayAssignments = new int[parentAssignments.size()];
                for(int i = 0; i < parentAssignments.size(); i++) arrayAssignments[i] = parentAssignments.get(i);
                parent.assignments = arrayAssignments;
                int parentPosIndex;
                if(level-1 > 0)
                    parentPosIndex = nodesByLevel.get(level-1).get(0).pos; // get the pos by looking at the first node
                else
                    parentPosIndex = -1;
                parent.pos = parentPosIndex;
                // don't bother actually updating the bounds, just add to nodes by level
                if(!nodesByLevel.containsKey(level-1))
                    nodesByLevel.put(level-1, new ArrayList<>());
                nodesByLevel.get(level - 1).add(parent);
            }
        }
        // the root node won't have any brethren
        return true;
    }

    public static boolean isPotentialChildOf(int[] parent, int[] potentialChild){
        for(int i =0; i < parent.length; i++){
            if (parent[i] != -1 && potentialChild[i] != parent[i]){
                return false;
            }
        }
        return true;
    }

    /*
    public static void recalculateScoresForNode(MultiSequenceSHARKStarNode node, RCs seqRCs, ObjectPool<MultiSequenceSHARKStarBound.ScoreContext> contexts, ConfAnalyzer confAnalyzer){
        try (ObjectPool.Checkout<MultiSequenceSHARKStarBound.ScoreContext> checkout = contexts.autoCheckout()) {
            MultiSequenceSHARKStarBound.ScoreContext context = checkout.get();
            node.index(context.index);
            double partialLB = context.partialConfLBScorer.calc(context.index, seqRCs);
            double partialUB = context.partialConfUBScorer.calc(context.index, seqRCs);
            double unassignLB = context.unassignedConfLBScorer.calc(context.index, seqRCs);
            double unassignUB = context.unassignedConfUBScorer.calc(context.index, seqRCs);
            ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, partialLB+unassignLB);
            ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
            double minimizedEnergy = analysis.epmol.energy;

            System.out.println(String.format("Recalculating %s: gscores [%.3f, %.3f], hscores [%.3f, %.3f], minE %.3f",
                    node.confToString(),
                    partialLB,
                    partialUB,
                    unassignLB,
                    unassignUB,
                    minimizedEnergy));
        }
    }

     */

    /*
    public static BigDecimal getZErrorReductionFromCorrections(SHARKStarQueue_refactor queue, Sequence seq){
        return queue.parallelStream()
                .map( (node) -> new BigMath(PartitionFunction.decimalPrecision).set(0.0)
                        .add(bc.calc(node.getPartialConfLB() + node.getUnassignedConfLB(seq)))
                        .sub(bc.calc(node.getPartialConfLB() + node.getUnassignedConfLB(seq) + node.getHOTCorrectionLB()))
                        .add(bc.calc(node.getPartialConfUB() + node.getUnassignedConfUB(seq) + node.getHOTCorrectionUB()))
                        .sub(bc.calc(node.getPartialConfUB() + node.getUnassignedConfUB(seq)))
                        .get()
                ).reduce(BigDecimal.ZERO, BigDecimal::add);
    }

    public static BigDecimal getZErrorReductionFromLowerCorrections(SHARKStarQueue_refactor queue, Sequence seq){
        return queue.parallelStream()
                .map( (node) -> new BigMath(PartitionFunction.decimalPrecision).set(0.0)
                        .add(bc.calc(node.getPartialConfLB() + node.getUnassignedConfLB(seq)))
                        .sub(bc.calc(node.getPartialConfLB() + node.getUnassignedConfLB(seq) + node.getHOTCorrectionLB()))
                        .get()
                ).reduce(BigDecimal.ZERO, BigDecimal::add);
    }

    public static BigDecimal getZErrorReductionFromUpperCorrections(SHARKStarQueue_refactor queue, Sequence seq){
        return queue.parallelStream()
                .map( (node) -> new BigMath(PartitionFunction.decimalPrecision).set(0.0)
                        .add(bc.calc(node.getPartialConfUB() + node.getUnassignedConfUB(seq) + node.getHOTCorrectionUB()))
                        .sub(bc.calc(node.getPartialConfUB() + node.getUnassignedConfUB(seq)))
                        .get()
                ).reduce(BigDecimal.ZERO, BigDecimal::add);
    }

    public static BigDecimal getZErrorReductionFromMinimizations(SHARKStarQueue_refactor queue, Sequence seq){
        return queue.parallelStream()
                .map( (node) -> {
                    if(node.isMinimized()){
                        return new BigMath(PartitionFunction.decimalPrecision).set(0.0)
                                .add(bc.calc(node.getPartialConfLB() + node.getUnassignedConfLB(seq) + node.getHOTCorrectionLB()))
                                .sub(bc.calc(node.getPartialConfUB() + node.getUnassignedConfUB(seq) + node.getHOTCorrectionUB()))
                                        .get();
                            }else{
                                return BigDecimal.ZERO;
                            }
                        }
                ).reduce(BigDecimal.ZERO, BigDecimal::add);
    }

    public static double getSumEnergyUpperCorrection(SHARKStarQueue_refactor queue, Sequence seq){
        return queue.parallelStream()
                .map(SHARKStarNode::getHOTCorrectionUB)
                .reduce(0.0, Double::sum);
    }

     */
}
