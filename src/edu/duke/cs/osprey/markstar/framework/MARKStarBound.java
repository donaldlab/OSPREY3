package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.*;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.astar.conf.RCs;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

public class MARKStarBound implements PartitionFunction {

    private double targetEpsilon = 1;
    private boolean debug = false;
    private Status status = null;
    private PartitionFunction.Values values = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;

    // the number of full conformations scored OR energied
    private int numConfsScored = 0;

    private boolean printMinimizedConfs;

    // Overwrite the computeUpperBound and computeLowerBound methods
    public static class Values extends PartitionFunction.Values {

        @Override
        public BigDecimal calcUpperBound() {
            return pstar;
        }

        @Override
        public BigDecimal calcLowerBound() {
            return qstar;
        }

        @Override
        public double getEffectiveEpsilon() {
            return MathTools.bigDivide(qstar, pstar, decimalPrecision).doubleValue();
        }
    }

    public void setReportProgress(boolean showPfuncProgress) {
        this.printMinimizedConfs = true;
    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    @Override
    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        status = Status.Estimating;
        values = Values.makeFullRange();
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
        status = Status.Estimating;
        values = Values.makeFullRange();
    }


    @Override
    public Status getStatus() {
        return status;
    }

    @Override
    public PartitionFunction.Values getValues() {
        return values;
    }

    @Override
    public int getParallelism() {
        return 0;
    }

    @Override
    public int getNumConfsEvaluated() {
        return numConfsEnergied;
    }

    @Override
    public void compute(int maxNumConfs) {
        compute();

    }

    private void debugPrint(String s) {
        if(debug)
            System.out.println(s);
    }

    public void compute() {

        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;
        while (epsilonBound > targetEpsilon) {
            debugPrint("Tightening from epsilon of "+epsilonBound);
            tightenBound();
            debugPrint("Errorbound is now "+epsilonBound);
            if(lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("Error. Bounds got looser.");
                //System.exit(-1);
            }
            lastEps = epsilonBound;
        }
        status = Status.Estimated;
        values.qstar = rootNode.getLowerBound();
        values.pstar = rootNode.getUpperBound();
        values.qprime= rootNode.getUpperBound();
        rootNode.printTree();
    }

    @Override
    public PartitionFunction.Result makeResult() {
        PartitionFunction.Result result = new Result(getStatus(), getValues(), getNumConfsEvaluated());
        result.setNumConfsLooked(numConfsScored);
        return result;
    }

    public static class Builder {

        /** The energy matrix to use for upper bound scoring. */
        private EnergyMatrix rigidEmat;
        /** The energy matrix to use for calculating minimized energies. */
        private EnergyMatrix minimizingEmat;
        /** The ConfEnergyCalculator to use for n-body minimizied energies */
        private ConfEnergyCalculator minimizingConfEcalc;

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

        public Builder(EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc, SimpleConfSpace confSpace) {
            this(rigidEmat, minimizingEmat, minimizingConfEcalc, new RCs(confSpace));
        }

        public Builder(EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc, PruningMatrix pmat) {
            this(rigidEmat, minimizingEmat, minimizingConfEcalc, new RCs(pmat));
        }

        public Builder(EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc, RCs rcs) {
            this.rigidEmat = rigidEmat;
            this.minimizingEmat = minimizingEmat;
            this.minimizingConfEcalc = minimizingConfEcalc;
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
            this.gscorer = new PairwiseGScorer(minimizingEmat);
            this.hscorer = new TraditionalPairwiseHScorer(minimizingEmat, rcs);
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

        public Builder setReportProgress(boolean val) {
            reportProgress = val;
            return this;
        }

        public Builder setPruner(AStarPruner val) {
            pruner = val;
            return this;
        }

        public MARKStarBound build() {
            ConfAStarTree tree = new ConfAStarTree.Builder(minimizingEmat, rcs).build();
            if (reportProgress) {
                tree.initProgress();
            }
            return new MARKStarBound(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc ,rcs);
        }
    }

    /**
     * TODO: 1. Make MARKStarBounds use and update a queue.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */
    // We keep track of the root node for computing our K* bounds
    private MARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    private final PriorityQueue<MARKStarNode> queue;
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

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat, ConfEnergyCalculator minimizingConfEcalc, RCs rcs) {
        this.queue = new PriorityQueue<>();
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        hscorerFactory = (emats) -> new TraditionalPairwiseHScorer(emats, rcs); //MPLPPairwiseHScorer(updater, emats, 50, 0.03);

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs, gscorerFactory, hscorerFactory, true);
        queue.add(rootNode);
        updateBound();
        confIndex = new ConfIndex(rcs.getNumPos());
        this.RCs = rcs;
        this.order = new DynamicHMeanAStarOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat),hscorerFactory.make(minimizingEmat));
        this.pruner = null;

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.gscorer = gscorerFactory.make(minimizingEmat);
            context.hscorer = hscorerFactory.make(minimizingEmat);
            context.rigidscorer = gscorerFactory.make(rigidEmat);
            /** These scoreres should match the scorers in the MARKStarNode root - they perform the same calculations**/
            context.negatedhscorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, rigidEmat)); //this is used for upper bounds, so we want it rigid
            //context.negatedhscorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, minimizingEmat));
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        setParallelism(null);
    }

    private static class ScoreContext {
        public ConfIndex index;
        public AStarScorer gscorer;
        public AStarScorer hscorer;
        public AStarScorer negatedhscorer;
        public AStarScorer rigidscorer;
        public ConfEnergyCalculator ecalc;
    }



    public void setParallelism(Parallelism val) {

        if (val == null) {
            val = Parallelism.makeCpu(1);
        }

        parallelism = val;
        tasks = parallelism.makeTaskExecutor(1000);
        contexts.allocate(parallelism.getParallelism());
    }

    private void debugEpsilon(double curEpsilon) {
        if(debug && curEpsilon < epsilonBound) {
            System.err.println("Epsilon just got bigger.");
        }
    }


    public void tightenBound(){
        debugPrint("Current overall error bound: "+epsilonBound);
        if(queue.isEmpty()) {
            debugPrint("Out of conformations.");
        }
        MARKStarNode curNode = queue.poll();
        Node node = curNode.getConfSearchNode();
        debugPrint("Processing Node: "+node.toString()) ;

        //If the child is a leaf, calculate n-body minimized energies
        if(node.getLevel() == RCs.getNumPos() && !node.isMinimized()){
            tasks.submit(() -> {
                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    ScoreContext context = checkout.get();
                    ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, -node.getConfLowerBound());
                    ConfSearch.EnergiedConf econf = context.ecalc.calcEnergy(conf);
                    numConfsEnergied++;
                    //Assign true energies to the subtreeLowerBound and subtreeUpperBound
                    double energy = econf.getEnergy();
                    curNode.setBoundsFromConfLowerAndUpper(econf.getEnergy(), econf.getEnergy());
                    node.gscore = econf.getEnergy();
                    if(true &&
                            (energy < -node.getMinScore() || energy > -node.getMaxScore())) {
                        System.err.println("Bounds are incorrect:" + (-node.getMinScore()) + "!< " + energy + " or " + energy
                                + " !<" + (-node.getMaxScore()) + " Aborting.");
                        System.exit(-1);
                    }
                    String out = "Energy = " + String.format("%6.3e", energy)+", ["+(node.getConfLowerBound())+","+(node.getConfUpperBound())+"]";
                    debugPrint(out);

                    if (printMinimizedConfs) {
                        System.out.println("["+SimpleConfSpace.formatConfRCs(node.assignments)+"]" +String.format("conf:%4d, score:%12.6f, energy:%12.6f",
                                numConfsEnergied, econf.getScore(), econf.getEnergy()
                        )
                        +", bounds: "+epsilonBound);
                    }
                }
                return null;
            },
            // Dummy function. We're not doing anything here.
            (Node child)->{});
            tasks.waitForFinish();
            updateBound();
            return;
        }

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

                    node.index(context.index);
                    Node child = node.assign(nextPos, nextRc);

                    // score the child node differentially against the parent node
                    if(child.getLevel() < RCs.getNumPos()) {
                        double diff = context.gscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        child.gscore = diff;
                        double confLowerBound = child.gscore + hdiff;
                        double confUpperbound = rigiddiff + maxhdiff;
                        child.computeNumConformations(RCs);
                        child.setBoundsFromConfLowerAndUpper(confLowerBound,confUpperbound);
                        if(debug)
                            System.out.println(child);
                    }
                    if(child.getLevel() == RCs.getNumPos()) {
                        double confPairwiseLower = context.gscorer.calc(context.index.assign(nextPos, nextRc), RCs);
                        double confRigid = context.rigidscorer.calc(context.index.assign(nextPos, nextRc), RCs);
                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        if(confPairwiseLower > confRigid) {
                            System.err.println("Our bounds are not tight. Lower bound is " + (confPairwiseLower - confRigid) + " higher");
                            double temp = confRigid;
                            confRigid = confPairwiseLower;
                            confPairwiseLower = temp;
                        }
                        child.setBoundsFromConfLowerAndUpper(confPairwiseLower,confRigid);
                        child.gscore = child.getConfLowerBound();

                        numConfsScored++;
                        if(debug)
                            System.out.println(child);
                    }



                    return child;
                }

            }, (Node child) -> {

                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                // collect the possible children
                if (child.getScore() < Double.POSITIVE_INFINITY) {
                    children.add(MARKStarNodeChild);
                }
                if(!child.isMinimized())
                    queue.add(MARKStarNodeChild);
                else
                    MARKStarNodeChild.computeEpsilonErrorBounds();
            });
        }
        tasks.waitForFinish();
        curNode.markUpdated();
        updateBound();
    }

    private void updateBound() {
        double curEpsilon = epsilonBound;
        epsilonBound = rootNode.computeEpsilonErrorBounds();
        debugEpsilon(curEpsilon);
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
