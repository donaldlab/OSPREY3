package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.*;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.astar.conf.RCs;

import java.util.*;

public class MARKStarBound {

    public static class Builder {

        /** The energy matrix to use for pairwise residue conformation energies. */
        private EnergyMatrix emat;

        private RCs rcs;
        private AStarOrder order = null;
        private AStarScorer gscorer = null;
        private AStarScorer hscorer = null;
        private boolean showProgress = false;
        private ConfAStarFactory factory = new LinkedConfAStarFactory();
        private SimpleConfSpace confSpace = null;
        private AStarPruner pruner = null;
        private ConfAStarTree.Builder aStarBuilder;
        private boolean reportProgress = false;

        public Builder(EnergyMatrix emat, SimpleConfSpace confSpace) {
            this(emat, new RCs(confSpace));
        }

        public Builder(EnergyMatrix emat, PruningMatrix pmat) {
            this(emat, new RCs(pmat));
        }

        public Builder(EnergyMatrix emat, RCs rcs) {
            this.emat = emat;
            this.rcs = rcs;

            // Jeff: MPLP is dramatically faster for large A* searches
            // and for small searches, who cares how fast A* is,
            // so I think it makes a good default for all cases
            setMPLP();
        }

        public Builder setCustom(AStarOrder order, AStarScorer gscorer, AStarScorer hscorer) {
            aStarBuilder.setCustom(order, gscorer, hscorer);
            this.order = order;
            this.gscorer = gscorer;
            this.hscorer = hscorer;
            return this;
        }

        /**
         * Uses the traditional estimation function to guide the tree search.
         * {@cite Leach1998 Leach, A.R. and Lemon, A.P., 1998. Exploring the conformational
         * space of protein side chains using dead-end elimination and the A* algorithm.
         * Proteins Structure Function and Genetics, 33(2), pp.227-239.}
         */
        public Builder setTraditional() {
            aStarBuilder.setTraditional();
            this.order = new DynamicHMeanAStarOrder();
            this.gscorer = new PairwiseGScorer(emat);
            this.hscorer = new TraditionalPairwiseHScorer(emat, rcs);
            return this;
        }

        /**
         * Creates an A* search using a newer estimation function based on Max Product Linear
         * Programming (MPLP).
         * {@cite Globerson2008 Globerson, A. and Jaakkola, T.S., 2008. Fixing max-product: Convergent message passing
         * algorithms for MAP LP-relaxations. In Advances in neural information processing systems (pp. 553-560).}
         *
         * For large designs, this A* implementation can be dramatically faster than the traditional
         * one, and often require much less memory too.
         */
        public Builder setMPLP() {
            aStarBuilder.setMPLP();
            return this;
        }

        public Builder setMPLP(ConfAStarTree.MPLPBuilder builder) {
            aStarBuilder.setMPLP(builder);
            return this;
        }

        /**
         * Use external memory (eg, disk, SSD, NAS) when large A* searches
         * cannot fit in internal memory (eg, RAM).
         *
         * Use {@link ExternalMemory#setInternalLimit} to set the amount of fixed internal memory
         * and {@link ExternalMemory#setTempDir} to set the file path for external memory.
         */
        public Builder useExternalMemory() {
            ExternalMemory.checkInternalLimitSet();
            factory = new EMConfAStarFactory();
            return this;
        }

        public Builder setShowProgress(boolean val) {
            showProgress = val;
            return this;
        }

        public Builder setPruner(AStarPruner val) {
            pruner = val;
            return this;
        }

        public MARKStarBound build() {
            ConfAStarTree tree = new ConfAStarTree.Builder(emat, rcs).build();
            if (showProgress) {
                tree.initProgress();
            }
            return new MARKStarBound(confSpace, emat, rcs);
        }
    }

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

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix emat, RCs rcs) {
        this.queue = Queue.PriorityFactory.of(null);

        rootNode = MARKStarNode.makeRoot(confSpace, emat, rcs, true);
        queue.push(rootNode);
        confIndex = new ConfIndex(rcs.getNumPos());
        this.RCs = rcs;
        this.order = null;
        this.pruner = null;
    }

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
