package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;

import java.util.*;

public class MARKStarBound {

    /**
     * TODO: 1. Make MARKStarBounds use and update a queue.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */
    // We keep track of the root node for computing our K* bounds
    private MARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    private final Queue<MARKStarNode> queue;
    private double epsilonBound = Double.POSITIVE_INFINITY;
    private boolean boundChanged = false;
    private ConfIndex confIndex;
    public final AStarOrder order;
    public final AStarPruner pruner;
    private RCs RCs;
    private Parallelism parallelism;
    private TaskExecutor tasks;
    private ObjectPool<ScoreContext> contexts;

    private static class ScoreContext {
        public ConfIndex index;
        public AStarScorer gscorer;
        public AStarScorer hscorer;
    }

    public MARKStarBound(SimpleConfSpace confSpace, AStarOrder order,
                         AStarPruner pruner, EnergyMatrix emat, RCs rcs, boolean reportProgress){
        this.queue = Queue.PriorityFactory.of(null);

        rootNode = MARKStarNode.makeRoot(confSpace, emat, rcs, reportProgress);
        queue.push(rootNode);
        confIndex = new ConfIndex(rcs.getNumPos());
        this.order = order;
        this.RCs = rcs;
        this.pruner = pruner;
    }

    public double getBound(){
        if(boundChanged)
            recomputeBound(rootNode);
        return epsilonBound;
    }

    private void recomputeBound(MARKStarNode rootNode) {
    }

    public void setParallelism(Parallelism val) {

        if (val == null) {
            val = Parallelism.makeCpu(1);
        }

        parallelism = val;
        tasks = parallelism.makeTaskExecutor(1000);
        contexts.allocate(parallelism.getParallelism());
    }
    public void tightenBound(){
        MARKStarNode MARKStarNode = queue.poll();
        ConfAStarNode node = MARKStarNode.getConfSearchNode();
        RCs rcs = RCs;

        // which pos to expand next?
        int numChildren = 0;
        node.index(confIndex);
        int nextPos = order.getNextPos(confIndex, RCs);
        assert (!confIndex.isDefined(nextPos));
        assert (confIndex.isUndefined(nextPos));

        // score child nodes with tasks (possibly in parallel)
        List<MARKStarNode> children = new ArrayList<>();
        for (int nextRc : rcs.get(nextPos)) {

            if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                continue;
            }

            // if this child was pruned dynamically, then don't score it
            if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                continue;
            }

            tasks.submit(() -> {

                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    ScoreContext context = checkout.get();

                    // score the child node differentially against the parent node
                    node.index(context.index);
                    ConfAStarNode child = node.assign(nextPos, nextRc);
                    child.setGScore(context.gscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
                    child.setHScore(context.hscorer.calcDifferential(context.index, rcs, nextPos, nextRc));
                    return child;
                }

            }, (ConfAStarNode child) -> {

                // collect the possible children
                if (child.getScore() < Double.POSITIVE_INFINITY) {
                    MARKStarNode MARKStarNodeChild = MARKStarNode.makeChild(child);
                    children.add(MARKStarNodeChild);
                }
            });
        }
        tasks.waitForFinish();
        numChildren += children.size();
        queue.pushAll(children);
    }

    private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {

        // do we even have pruned pairs?
        PruningMatrix pmat = RCs.getPruneMat();
        if (pmat == null) {
            return false;
        }

        for (int i = 0; i < confIndex.numDefined; i++) {
            int pos = confIndex.definedPos[i];
            int rc = confIndex.definedRCs[i];
            assert (pos != nextPos || rc != nextRc);
            if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
                return true;
            }
        }
        return false;
    }
}
