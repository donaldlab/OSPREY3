package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.*;
import edu.duke.cs.osprey.astar.conf.linked.LinkedConfAStarFactory;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.order.UpperLowerAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.*;
import edu.duke.cs.osprey.externalMemory.EMConfAStarFactory;
import edu.duke.cs.osprey.externalMemory.ExternalMemory;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.tools.Stopwatch;
import javafx.util.Pair;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;

public class MARKStarBound implements PartitionFunction {

    private double targetEpsilon = 1;
    private boolean debug = false;
    private Status status = null;
    private MARKStarBound.Values values = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;
    // max confs minimized, -1 means infinite.
    private int maxNumConfs = -1;

    // the number of full conformations scored OR energied
    private int numConfsScored = 0;

    private boolean printMinimizedConfs;
    private MARKStarProgress progress;
    public String stateName = String.format("%4f",Math.random());
    private int numPartialMinimizations;

    // Overwrite the computeUpperBound and computeLowerBound methods
    public static class Values extends PartitionFunction.Values {

        public Values ()
        {
            pstar = MathTools.BigPositiveInfinity;
        }
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
            return MathTools.bigDivide(pstar.subtract(qstar), pstar, decimalPrecision).doubleValue();
        }
    }

    public void setReportProgress(boolean showPfuncProgress) {
        this.printMinimizedConfs = true;
    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    public void setMaxNumConfs(int maxNumConfs) {
        this.maxNumConfs = maxNumConfs;
    }

    @Override
    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        status = Status.Estimating;
        values = new MARKStarBound.Values();
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
        status = Status.Estimating;
        values = new MARKStarBound.Values();
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
        return numConfsEnergied + numPartialMinimizations;
    }

    @Override
    public int getNumConfsScored() {
        return numConfsScored;
    }

    @Override
    public void compute(int maxNumConfs) {
        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;
        while (epsilonBound > targetEpsilon &&
                (maxNumConfs < 0 || numConfsEnergied < maxNumConfs)) {
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
        rootNode.printTree(stateName);
    }

    private void debugPrint(String s) {
        if(debug)
            System.out.println(s);
    }

    public void compute() {
        compute(-1);
    }

    @Override
    public PartitionFunction.Result makeResult() {
        PartitionFunction.Result result = new PartitionFunction.Result(getStatus(), getValues(), getNumConfsEvaluated(), numConfsScored, rootNode.getNumConfs());
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
        private Parallelism parallelism;
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

        public Builder setParalllelism(Parallelism val) {
            parallelism = val;
            return this;
        }
        public MARKStarBound build() {
            ConfAStarTree tree = new ConfAStarTree.Builder(minimizingEmat, rcs).build();
            if (reportProgress) {
                tree.initProgress();
            }
            return new MARKStarBound(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc ,rcs, parallelism);
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
    // Heap of nodes (pointer copies) of upcoming minimization
    private final Queue<MARKStarNode> minimizationQueue;
    private double epsilonBound = Double.POSITIVE_INFINITY;
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
    private int stepSize;
    public static final int MAX_STEP_SIZE = 20;

    public static final int MAX_CONFSPACE_FRACTION = 1000000;
    public static final double MINIMIZATION_FACTOR = 0.1;
    public boolean reduceMinimizations = false;
    private ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc = null;
    private Stopwatch stopwatch = new Stopwatch().start();
    double msRunning = 0;
    private static final int ReportIntervalMs = 10 * 1000; // TODO: make configurable

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.queue = new PriorityQueue<>();
        this.minimizationQueue = new PriorityQueue<>();
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        hscorerFactory = (emats) -> new MPLPPairwiseHScorer(updater, emats, 50, 0.03);
        this.correctionMatrix = new UpdatingEnergyMatrix(confSpace, minimizingEmat);

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs, gscorerFactory, hscorerFactory, true);
        queue.add(rootNode);
        updateBound();
        confIndex = new ConfIndex(rcs.getNumPos());
        this.minimizingEmat = minimizingEmat;
        this.RCs = rcs;
        this.order = new UpperLowerAStarOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat),hscorerFactory.make(minimizingEmat));
        this.pruner = null;
        stepSize = Math.min(MAX_STEP_SIZE,
                Math.max(1,RCs.getNumConformations().divide(new BigInteger(""+MAX_CONFSPACE_FRACTION)).intValue()));

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.gscorer = gscorerFactory.make(minimizingEmat);
            context.hscorer = hscorerFactory.make(minimizingEmat);
            context.rigidscorer = gscorerFactory.make(rigidEmat);
            /** These scoreres should match the scorers in the MARKStarNode root - they perform the same calculations**/
            context.negatedhscorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, rigidEmat)); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        progress = new MARKStarProgress(RCs.getNumPos());
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        setParallelism(parallelism);
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


    public void tightenBound() {
        debugPrint("Current overall error bound: "+epsilonBound);
        if(queue.isEmpty()) {
            debugPrint("Out of conformations.");
        }
        int numStepsThisLoop = 0;
        boolean minimizing = false;
        List<MARKStarNode> newNodes = new ArrayList<>();
        List<MARKStarNode> newNodesToMinimize = new ArrayList<>();
        while(!queue.isEmpty() && epsilonBound > targetEpsilon) {
                if(epsilonBound <= targetEpsilon)
                    break;
                if(numStepsThisLoop >= stepSize || (minimizing && reduceMinimizations))
                    break;
                synchronized (this)
                {
                    numStepsThisLoop++;
                }
                MARKStarNode curNode = queue.poll();
                Node node = curNode.getConfSearchNode();
                debugPrint("Processing Node: " + node.toString());
                if(numStepsThisLoop < 2 &&!reduceMinimizations
                        && curNode.getConfSearchNode().getSubtreeUpperBound().compareTo(new BigDecimal(100))<1)
                {
                    System.err.println("Node error is insignificant. Why is this happening? Aren't we done?");
                }

                //If the child is a leaf, calculate n-body minimized energies
                if (node.getLevel() == RCs.getNumPos() && !node.isMinimized()) {
                    minimizing = true;
                    tasks.submit(() -> {
                                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                                    ScoreContext context = checkout.get();
                                    node.index(context.index);
                                    double confCorrection = correctionMatrix.confE(node.assignments);
                                    if(node.getConfLowerBound()!= confCorrection) {
                                        debugPrint("Correcting :["+SimpleConfSpace.formatConfRCs(node.assignments)
                                                +":"+node.gscore+"] down to "+confCorrection);
                                        node.gscore = confCorrection;
                                        if(confCorrection > node.rigidScore)
                                            System.err.println("Overcorrected: "+confCorrection + " > "+node.rigidScore);
                                        node.setBoundsFromConfLowerAndUpper(confCorrection, node.rigidScore);
                                        curNode.markUpdated();
                                        updateBound();
                                        return null;
                                    }
                                    ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound());
                                    ConfSearch.EnergiedConf econf = context.ecalc.calcEnergy(conf);
                                    double energy = econf.getEnergy();
                                    double newConfUpper = econf.getEnergy();
                                    double newConfLower = econf.getEnergy();
                                    if(energy < node.getConfLowerBound()) {
                                        System.err.println("Bounds are incorrect:" + (node.getConfLowerBound()) + " > "
                                                + energy);
                                        if (energy < 10)
                                            System.err.println("The bounds are probably wrong.");
                                            //System.exit(-1);
                                    }
                                    if (energy > node.getConfUpperBound()) {
                                        System.err.println("Upper bounds got worse after minimization:" + energy
                                                + " > " + (node.getConfUpperBound())+". Rejecting minimized energy.");
                                        System.err.println("Node info: "+node);
                                        double minimized = context.ecalc.calcEnergy(conf).getEnergy();
                                        ConfIndex nodeIndex = new ConfIndex(node.assignments.length);
                                        node.index(nodeIndex);
                                        double rigid = context.rigidscorer.calc(nodeIndex, RCs);
                                        double pairiwseMin = context.gscorer.calc(nodeIndex, RCs);

                                        newConfUpper = node.getConfUpperBound();
                                        newConfLower = node.getConfUpperBound();
                                    }
                                    computeEnergyCorrection(conf, context.gscorer, context.ecalc);
                                    curNode.setBoundsFromConfLowerAndUpper(newConfLower,newConfUpper);
                                    node.gscore = newConfLower;
                                    String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound()) + "," + (node.getConfUpperBound()) + "]";
                                    debugPrint(out);
                                    updateBound();
                                    synchronized(this) {
                                        numConfsEnergied++;
                                    }
                                    //Assign true energies to the subtreeLowerBound and subtreeUpperBound
                                    if (printMinimizedConfs) {
                                        System.out.println("[" + SimpleConfSpace.formatConfRCs(node.assignments) + "]" + String.format("conf:%4d, score:%12.6f, energy:%12.6f",
                                                numConfsEnergied, econf.getScore(), newConfLower
                                        )
                                                + ", bounds: " + epsilonBound);
                                    }
                                }
                                return null;
                            },
                            // Dummy function. We're not doing anything here.
                            (Node child) -> {
                                progress.reportLeafNode(node.gscore, queue.size(), epsilonBound);
                                if(!node.isMinimized())
                                    newNodes.add(curNode);

                            });
                    if(epsilonBound <= targetEpsilon)
                        return;
                    continue;
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
                            if (child.getLevel() < RCs.getNumPos()) {
                                double diff = context.gscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                child.gscore = diff;
                                //Correct for incorrect gscore.
                                rigiddiff=rigiddiff-node.gscore+node.rigidScore;
                                child.rigidScore = rigiddiff;

                                double confLowerBound = child.gscore + hdiff;
                                double confUpperbound = rigiddiff + maxhdiff;
                                child.computeNumConformations(RCs);
                                child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                                if (debug)
                                    System.out.println(child);
                                progress.reportInternalNode(child.level,child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                            }
                            if (child.getLevel() == RCs.getNumPos()) {
                                double confPairwiseLower = context.gscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                                confRigid=confRigid-node.gscore+node.rigidScore;

                                child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                                if (confPairwiseLower > confRigid) {
                                    System.err.println("Our bounds are not tight. Lower bound is " + (confPairwiseLower - confRigid) + " higher");
                                    double temp = confRigid;
                                    confRigid = confPairwiseLower;
                                    confPairwiseLower = temp;
                                }
                                double confCorrection = correctionMatrix.confE(child.assignments);
                                double lowerbound = minimizingEmat.confE(child.assignments);
                                if(lowerbound != confCorrection)
                                    debugPrint("Correcting node "+SimpleConfSpace.formatConfRCs(child.assignments)
                                    +":"+lowerbound+"->"+confCorrection);
                                checkBounds(confCorrection,confRigid);
                                correctionMatrix.confE(child.assignments);
                                child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                                child.gscore = child.getConfLowerBound();
                                child.rigidScore = confRigid;
                                if(reduceMinimizations)
                                    child.setMinimizationRatio(MINIMIZATION_FACTOR/(RCs.getNumPos()*RCs.getNumPos()));
                                numConfsScored++;
                                if (debug)
                                    System.out.println(child);
                                progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                            }


                            return child;
                        }

                    }, (Node child) -> {
                        if(Double.isNaN(child.rigidScore))
                            System.out.println("Huh!?");
                        MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                        synchronized (this) {
                            // collect the possible children
                            if (child.getScore() < Double.POSITIVE_INFINITY) {
                                children.add(MARKStarNodeChild);
                            }
                            if (!child.isMinimized()) {
                                newNodes.add(MARKStarNodeChild);
                                newNodesToMinimize.add(MARKStarNodeChild);
                            }
                            else
                                MARKStarNodeChild.computeEpsilonErrorBounds();

                        }
                        curNode.markUpdated();
                    });
                }
        }

        tasks.waitForFinish();
        queue.addAll(newNodes);
        minimizationQueue.addAll(newNodesToMinimize);
        if(!reduceMinimizations || minimizationQueue.size() > 200)
            processPreminimization(minimizingEcalc);
        debugHeap();
        updateBound();


    }


    private void checkBounds(double lower, double upper)
    {
        if (lower > upper)
            System.err.println("Bounds incorrect.");
    }

    private void computeEnergyCorrection(ConfSearch.ScoredConf conf, AStarScorer gscorer, ConfEnergyCalculator ecalc) {
        // TODO: Replace the sortedPairwiseTerms with an ArrayList<TupE>.
        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
        //System.out.println("Analysis:"+analysis);
        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition();
        //System.out.println("Energy Analysis: "+energyAnalysis);
        //System.out.println("Score Analysis: "+scoreAnalysis);
        EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);
        //System.out.println("Difference Analysis " + diff);
        List<Pair<Pair<Integer, Integer>, Double>> sortedPairwiseTerms = new ArrayList<>();
        for (int pos = 0; pos < diff.getNumPos(); pos++)
        {
            for (int rc = 0; rc < diff.getNumConfAtPos(pos); rc++)
            {
                for (int pos2 = 0; pos2 < diff.getNumPos(); pos2++)
                {
                    for (int rc2 = 0; rc2 < diff.getNumConfAtPos(pos2); rc2++)
                    {
                        if(pos >= pos2)
                            continue;
                        double sum = 0;
                        sum+=diff.getOneBody(pos, rc);
                        sum+=diff.getPairwise(pos, rc, pos2, rc2);
                        sum+=diff.getOneBody(pos2,rc2);
                        sortedPairwiseTerms.add(new Pair(new Pair(pos,pos2), sum));
                    }
                }
            }
        }
        Collections.sort(sortedPairwiseTerms, (a,b)->-Double.compare(a.getValue(),b.getValue()));

        //Collections.sort(sortedPairwiseTerms, Comparator.comparingDouble(Pair::getValue));
        double threshhold = 0.3;
        for(int i = 0; i < sortedPairwiseTerms.size(); i++)
        {
            Pair<Pair<Integer, Integer>, Double> pairEnergy = sortedPairwiseTerms.get(i);
            if(pairEnergy.getValue() < threshhold)
                continue;
            int pos1 = pairEnergy.getKey().getKey();
            int pos2 = pairEnergy.getKey().getValue();
            double maxCorrection = 0;
            int localMinimizations = 0;
            for(int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double correction = computeDifference(tuple, ecalc);
                localMinimizations++;
                debugPrint("Correction for "+tuple.stringListing()+":"+correction);
                if(correction > 0 )
                    correctionMatrix.setHigherOrder(tuple, correction);
                else
                    System.err.println("Positive correction for "+tuple.stringListing());
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        /* Starting from the largest-difference pairs, create triples and quads to find
         * tuples which correct the energy of the pair.
         */
    }




    private double computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        double tripleEnergy = ecalc.calcEnergy(tuple).energy;
        double lowerbound = minimizingEmat.getInternalEnergy(tuple);
        return tripleEnergy-lowerbound;
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for(int pos: positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    private void processPreminimization(ConfEnergyCalculator ecalc)
    {
        int maxMinimizations = 10;
        List<MARKStarNode> topConfs = getTopConfs(maxMinimizations);
        // Need at least two confs to do any partial preminimization
        if(topConfs.size() < 2)
            return;
        RCTuple lowestBoundTuple= topConfs.get(0).toTuple();
        RCTuple overlap = findLargestOverlap(lowestBoundTuple, topConfs, 4);
        //Only continue if we have something to minimize
        if(overlap.size() < 4 || correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        System.out.println("Computing correction for "+overlap.stringListing()+" penalty of "+(partiallyMinimizedLower-pairwiseLower));
        progress.reportPartialMinimization(1, epsilonBound);
        correctionMatrix.setHigherOrder(overlap, partiallyMinimizedLower-pairwiseLower);
        minimizationQueue.addAll(topConfs);
    }

    private List<MARKStarNode> getTopConfs(int numConfs) {
        List<MARKStarNode> topConfs = new ArrayList<>();
        while (topConfs.size() < numConfs&& !minimizationQueue.isEmpty()) {
            MARKStarNode nextLowestConf = minimizationQueue.poll();
            topConfs.add(nextLowestConf);
        }
        return topConfs;
    }


    private RCTuple findLargestOverlap(RCTuple conf, List<MARKStarNode> otherConfs, int minResidues) {
        RCTuple overlap = conf;
        for(MARKStarNode other: otherConfs) {
            overlap = overlap.intersect(other.toTuple());
            if(overlap.size() < minResidues)
                break;
        }
        return overlap;

    }

    private void debugHeap() {
        if(debug) {
            if(queue.isEmpty())
                return;
            debugPrint("A* Heap:");
            PriorityQueue<MARKStarNode> copy = new PriorityQueue<>();
            BigDecimal peek = queue.peek().getErrorBound();
            int numConfs = 0;
            while(!queue.isEmpty() && numConfs < 10) {
                MARKStarNode node = queue.poll();
                if(node.getConfSearchNode().isLeaf() || node.getErrorBound().multiply(new BigDecimal(1000)).compareTo(peek) >= 0)
                debugPrint(node.getConfSearchNode().toString() + "," + String.format("%12.6e", node.getErrorBound()));
                copy.add(node);
                numConfs++;
            }
            queue.addAll(copy);
        }
    }

    private void updateBound() {
        double curEpsilon = epsilonBound;
        synchronized (this) {
            epsilonBound = rootNode.computeEpsilonErrorBounds();
        }
        debugEpsilon(curEpsilon);
        //System.out.println("Current epsilon:"+epsilonBound);
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
