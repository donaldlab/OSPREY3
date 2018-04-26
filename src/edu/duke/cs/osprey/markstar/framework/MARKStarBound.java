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
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.astar.conf.RCs;

import java.math.BigDecimal;
import java.util.*;

public class MARKStarBound implements PartitionFunction {

    private double targetEpsilon = 1;

    public void setReportProgress(boolean showPfuncProgress) {
    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    @Override
    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
    }

    @Override
    public Status getStatus() {
        return null;
    }

    @Override
    public Values getValues() {
        return null;
    }

    @Override
    public int getParallelism() {
        return 0;
    }

    @Override
    public int getNumConfsEvaluated() {
        return 0;
    }

    @Override
    public void compute(int maxNumConfs) {

    }

    public void compute() {
        while (epsilonBound > targetEpsilon) {
            System.out.println("Tightening from epsilon of "+epsilonBound);
            tightenBound();
            System.out.println("Errorbound is now "+epsilonBound);
        }
    }

    public PartitionFunction.Result makeResult() {
        return null;
    }

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
    // TODO: Implement new AStarPruner for MARK*?
    public final AStarPruner pruner;
    private RCs RCs;
    private Parallelism parallelism;
    private TaskExecutor tasks;
    private ObjectPool<ScoreContext> contexts;
    private MARKStarNode.ScorerFactory gscorerFactory;
    private MARKStarNode.ScorerFactory hscorerFactory;

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix emat, RCs rcs) {
        this.queue = Queue.PriorityFactory.of(null);
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        hscorerFactory = (emats) -> new TraditionalPairwiseHScorer(emats, rcs); //MPLPPairwiseHScorer(updater, emats, 50, 0.03);

        rootNode = MARKStarNode.makeRoot(confSpace, emat, rcs, gscorerFactory, hscorerFactory, true);
        queue.push(rootNode);
        updateBound();
        confIndex = new ConfIndex(rcs.getNumPos());
        this.RCs = rcs;
        this.order = new DynamicHMeanAStarOrder();
        order.setScorers(gscorerFactory.make(emat),hscorerFactory.make(emat));
        this.pruner = null;

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.gscorer = gscorerFactory.make(emat);
            context.hscorer = hscorerFactory.make(emat);
            context.negatedhscorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, emat));
            return context;
        });

        setParallelism(null);
    }

    private static class ScoreContext {
        public ConfIndex index;
        public AStarScorer gscorer;
        public AStarScorer hscorer;
        public AStarScorer negatedhscorer;
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
        System.out.println("Current overall error bound: "+epsilonBound);
        if(queue.isEmpty()) {
            System.out.println("Out of conformations.");
        }
        MARKStarNode curNode = queue.poll();
        Node node = curNode.getConfSearchNode();
        System.out.println("Processing Node: "+node.toString()) ;

        // which pos to expand next?
        int numChildren = 0;
        node.index(confIndex);
        int nextPos = order.getNextPos(confIndex, RCs);
        assert (!confIndex.isDefined(nextPos));
        assert (confIndex.isUndefined(nextPos));

        // score child nodes with tasks (possibly in parallel)
        List<MARKStarNode> children = new ArrayList<>();
        for (int nextRc : RCs.get(nextPos)) {

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
                    Node child = node.assign(nextPos, nextRc);
                    //TODO: Change this code to do the right thing.
                    System.out.println("Expanded node:"+child.confToString());
                    double diff = context.gscorer.calcDifferential(context.index,RCs,nextPos, nextRc);
                    double hdiff = context.hscorer.calcDifferential(context.index,RCs,nextPos,nextRc);
                    double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    child.gscore = node.gscore + diff;
                    child.minHScore = child.gscore + hdiff;
                    child.maxHScore = child.gscore + maxhdiff;
                    System.out.println("g score:"+child.gscore+", min:"+child.minHScore+", max:"+child.maxHScore);
                    return child;
                }

            }, (Node child) -> {

                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                // collect the possible children
                if (child.getScore() < Double.POSITIVE_INFINITY) {
                    children.add(MARKStarNodeChild);
                }
            });
        }
        tasks.waitForFinish();
        updateBound();
        numChildren += children.size();
        for(MARKStarNode child: children) {
            //Only add partial conformations to the queue.
            if(child.level < RCs.getNumPos())
                queue.push(child);
        }
    }

    private void updateBound() {
        epsilonBound = rootNode.getErrorBound();
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
