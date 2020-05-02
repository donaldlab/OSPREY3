package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.tools.ObjectPool;
import org.jetbrains.annotations.NotNull;

import java.util.AbstractMap.SimpleEntry;
import java.util.List;
import java.util.stream.Collectors;

public class SHARKStarQueueDebugger {
    /**
     * A class to debug a SHARKStarQueue_refactor object
     */


    public SHARKStarQueueDebugger(){

    }

    /**
     * Checks for duplicates in the SHARKStarQueue
     *
     * @param queue The queue to check
     * @return true if there are duplicates, false otherwise
     */
    public static boolean hasDuplicateNodes(SHARKStarQueue_refactor queue){
        return queue.parallelStream().map( node -> new SimpleEntry<>(node, 1)) // map nodes to nodes with freq counter
                .collect(Collectors.groupingBy(e -> e.getKey().getAssignments(), Collectors.counting())) //collect into map with freqs
                .values().stream().map(e -> e>1).reduce(false, (a,b) -> a || b); // check whether any frequencies are greater than 1
    }

    /**
     * Prints the number of nodes at each level in queue
     *
     * @param queue The queue to check
     */
    public static void printLevelBreakdown(SHARKStarQueue_refactor queue){
        queue.parallelStream().map(node -> new SimpleEntry<>(node.getLevel(), 1)) // map nodes to level and freq counter
            .collect(Collectors.groupingBy(SimpleEntry::getKey, Collectors.counting()))
            .forEach((k, v) -> System.out.println(String.format("Level %d: %d nodes", k, v)));
    }

    /**
     * Prints the node information for sequence for each node at level
     *
     * @param queue The queue to check
     * @param level The level at which to print nodes
     * @param seq   The sequence for which to print information
     */
    public static void printNodesAtLevel(SHARKStarQueue_refactor queue, int level, Sequence seq){
        System.out.println(String.format("Printing all nodes at level %d:", level));
        getNodesAtLevel(queue, level).forEach(node -> System.out.println(String.format("%s: [%.3f, %.3f], correct: %.3f, minE: %.3f",
                        node.confToString(),
                        node.getFreeEnergyLB(seq),
                        node.getFreeEnergyUB(seq),
                        node.getHOTCorrection(),
                        node.getMinE()
                        )));
    }

    /**
     * Returns the nodes at level
     *
     * @param queue The queue to check
     * @param level The level at which to print nodes
     */
    public static List<SHARKStarNode> getNodesAtLevel(SHARKStarQueue_refactor queue, int level) {
        return queue.parallelStream().filter(node -> node.getLevel() == level)
                .collect(Collectors.toList());
    }

    public static void recalculateScoresForNode(SHARKStarNode node, RCs seqRCs, ObjectPool<MultiSequenceSHARKStarBound_refactor.ScoreContext> contexts){
        try (ObjectPool.Checkout<MultiSequenceSHARKStarBound_refactor.ScoreContext> checkout = contexts.autoCheckout()) {
            MultiSequenceSHARKStarBound_refactor.ScoreContext context = checkout.get();
            node.index(context.index);
            double partialLB = context.partialConfLBScorer.calc(context.index, seqRCs);
            double partialUB = context.partialConfUBScorer.calc(context.index, seqRCs);
            double unassignLB = context.unassignedConfLBScorer.calc(context.index, seqRCs);
            double unassignUB = context.unassignedConfUBScorer.calc(context.index, seqRCs);

            System.out.println(String.format("Recalculating %s: gscores [%.3f, %.3f], hscores [%.3f, %.3f]",
                    node.confToString(),
                    partialLB,
                    partialUB,
                    unassignLB,
                    unassignUB));
        }
    }
}
