/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
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
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupE;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.PriorityBlockingQueue;

public class MARKStarBoundAsync implements PartitionFunction {

    private double targetEpsilon = 1;
    private BigDecimal targetEpsilonAsBD = new BigDecimal(targetEpsilon);
    public boolean debug = false;
    public boolean profileOutput = false;
    private Status status = null;
    private Values values = null;

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
    private ArrayList<Integer> minList;
    private double internalTimeAverage;
    private double leafTimeAverage;
    private double cleanupTime;
    private boolean nonZeroLower;
    private TaskExecutor loopTasks;
    private TaskExecutor correctionTasks;

    public void setCorrections(UpdatingEnergyMatrix cachedCorrections) {
        correctionMatrix = cachedCorrections;
    }

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

    @Override
    public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon){
        state = new State(numConfsBeforePruning);
        init(targetEpsilon);
    }

    public void setRCs(RCs rcs) {
        RCs = rcs;
    }

    public void setReportProgress(boolean showPfuncProgress) {
        this.printMinimizedConfs = true;
    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    @Override
    public void setStabilityThreshold(BigDecimal threshhold) {
        stabilityThreshold = threshhold;
    }

    public void setMaxNumConfs(int maxNumConfs) {
        this.maxNumConfs = maxNumConfs;
    }

    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        this.targetEpsilonAsBD = new BigDecimal(targetEpsilon);
        status = Status.Estimating;
        values = new Values();
        state = new State(RCs.getNumConformations());
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
        status = Status.Estimating;
        values = new Values();
        state = new State(RCs.getNumConformations());
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

    public int getNumConfsScored() {
        return numConfsScored;
    }

    @Override
    public void compute(int maxNumConfs) {
        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;

        int previousConfCount = numConfsEnergied + numConfsScored + numPartialMinimizations;

        if(!nonZeroLower) {
            runUntilNonZero();
            updateBound();
        }
        while (epsilonBound > targetEpsilon &&
                (maxNumConfs < 0 || numConfsEnergied + numConfsScored + numPartialMinimizations - previousConfCount < maxNumConfs)) {
            debugPrint("Tightening from epsilon of "+epsilonBound);
            tightenBoundAsyncDirect();
            debugPrint("Errorbound is now "+epsilonBound);
            if(lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("Error. Bounds got looser.");
                //System.exit(-1);
            }
            lastEps = epsilonBound;
        }
        if(epsilonBound > targetEpsilon) {
            loopTasks.waitForFinish();
            minimizingEcalc.tasks.waitForFinish();
        }
        BigDecimal averageReduction = BigDecimal.ZERO;
        int totalMinimizations = numConfsEnergied + numPartialMinimizations;
        if(totalMinimizations> 0)
            averageReduction = cumulativeZCorrection
                .divide(new BigDecimal(totalMinimizations), new MathContext(BigDecimal.ROUND_HALF_UP));
        debugPrint(String.format("Average Z reduction per minimization: %12.6e",averageReduction));
        if(epsilonBound < targetEpsilon)
            status = Status.Estimated;
        values.qstar = rootNode.getLowerBound();
        values.pstar = rootNode.getUpperBound();
        values.qprime= rootNode.getUpperBound();
        //rootNode.printTree(stateName);
    }

    private void debugPrint(String s) {
        if(debug)
            System.out.println(s);
    }

    private void profilePrint(String s) {
        if(profileOutput)
            System.out.println(s);
    }

    public void compute() {
        compute(Integer.MAX_VALUE);
    }

    @Override
    public Result makeResult() {
        Result result = new Result(getStatus(), getValues(), getNumConfsEvaluated());
        return result;
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
    private ConfIndex confIndex;
    public final AStarOrder order;
    // TODO: Implement new AStarPruner for MARK*?
    public final AStarPruner pruner;
    private RCs RCs;
    private Parallelism parallelism;
    private TaskExecutor internalTasks;
    private TaskExecutor leafTasks;
    private TaskExecutor drillTasks;
    private ObjectPool<ScoreContext> contexts;
    private MARKStarNode.ScorerFactory gscorerFactory;
    private MARKStarNode.ScorerFactory hscorerFactory;

    public static final int MAX_CONFSPACE_FRACTION = 1000000;
    public static final double MINIMIZATION_FACTOR = 0.1;
    public boolean reduceMinimizations = true;
    private ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix rigidEmat;
    UpdatingEnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc = null;
    private Stopwatch stopwatch = new Stopwatch().start();
    BigDecimal cumulativeZCorrection = BigDecimal.ZERO;
    BigDecimal ZReductionFromMin = BigDecimal.ZERO;
    BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private boolean computedCorrections = false;
    private long loopPartialTime = 0;
    private Set<String> correctedTuples = Collections.synchronizedSet(new HashSet<String>());
    private State state;
    private BigDecimal stabilityThreshold;
    BlockingQueue<MARKStarNode> asyncQueue = new PriorityBlockingQueue<>();
    private double leafTimeSum = 0;
    private double internalTimeSum = 0;
    private int numLeavesScored = 0;
    private int numInternalScored = 0;
    private double bias = 1;

    public static MARKStarBoundFastQueues makeFromConfSpaceInfo(BBKStar.ConfSpaceInfo info, RCs rcs) {
        throw new UnsupportedOperationException("MARK* is not yet integrated into BBK*. Coming soon!");
        /*
        ConfEnergyCalculator minimizingConfEcalc = info.confEcalcMinimized;
        return new MARKStarBound(info.confSpace, info.ematRigid, info.ematMinimized, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
        */
    }

    public MARKStarBoundAsync(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                              ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.queue = new PriorityQueue<>();
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        hscorerFactory = (emats) -> new MPLPPairwiseHScorer(updater, emats, 50, 0.03);

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
                gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat),
                gscorerFactory.make(rigidEmat),
                new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);
        confIndex = new ConfIndex(rcs.getNumPos());
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.RCs = rcs;
        this.order = new UpperLowerAStarOrder();
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
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        progress = new MARKStarProgress(RCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        setParallelism(parallelism);
        updateBound();
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(),0));
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
        loopTasks = parallelism.makeTaskExecutor(10000);
        correctionTasks = parallelism.makeTaskExecutor(10000);
        contexts.allocate(parallelism.getParallelism());
    }

    private void debugEpsilon(double curEpsilon) {
        if(debug && curEpsilon < epsilonBound) {
            System.err.println("Epsilon just got bigger.");
        }
    }

    private boolean shouldMinimize(Node node) {
        return node.getLevel() == RCs.getNumPos() && !node.isMinimized();
    }

    private void recordCorrection(double lowerBound, double correction) {
        BigDecimal upper = bc.calc(lowerBound);
        BigDecimal corrected = bc.calc(lowerBound + correction);
        cumulativeZCorrection = cumulativeZCorrection.add(upper.subtract(corrected));
    }
    private void recordReduction(double score, double energy) {
        BigDecimal scoreWeight = bc.calc(score);
        BigDecimal energyWeight = bc.calc(energy);
        ZReductionFromMin = ZReductionFromMin.add(scoreWeight.subtract(energyWeight));

    }

    // We want to process internal nodes without worrying about the bound too much until we have
    // a nonzero lower bound. We have to have a nonzero lower bound, so we have to have at least
    // one node with a negative conf upper bound.
    private void runUntilNonZero() {
        System.out.println("Running until leaf is found...");
        double bestConfUpper = Double.POSITIVE_INFINITY;

        List<MARKStarNode> newNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        boundLowestBoundConfUnderNode(rootNode,newNodes);
        queue.addAll(newNodes);
        asyncQueue.addAll(newNodes);


        newNodes.clear();
        System.out.println("Found a leaf!");
        nonZeroLower = true;
    }


    private void tightenBoundAsyncDirect() {
        System.out.println(String.format("Current overall error bound: %12.10f, spread of [%12.6e, %12.6e]",epsilonBound, rootNode.getLowerBound(), rootNode.getUpperBound()));
        List<MARKStarNode> internalNodes = state.internalNodes;
        List<MARKStarNode> leafNodes = state.leafNodes;
        List<MARKStarNode> newNodes = new ArrayList<>();
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        int numNodes = 0;
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        while(!asyncQueue.isEmpty() && !asyncQueue.peek().getConfSearchNode().isLeaf()) {
            MARKStarNode internalNode = asyncQueue.poll();{
                Stopwatch internalTime = new Stopwatch();
                if(!MathTools.isGreaterThan(internalNode.getLowerBound(),BigDecimal.ZERO) &&
                    MathTools.isGreaterThan(
                            MathTools.bigDivide(internalNode.getUpperBound(),rootNode.getUpperBound(),
                                    PartitionFunction.decimalPrecision),
                            new BigDecimal(1-targetEpsilon))
                ) {
                    loopTasks.submit(() -> {
                        internalTime.start();
                        boundLowestBoundConfUnderNodeAsync(internalNode);
                        internalTime.stop();
                        return internalTime.getTimeS();
                    }, (timeS) -> {
                        internalNode.markUpdated();
                        asyncQueue.addAll(newNodes);
                        onFinishedInternal(timeS);
                    });
                }
                else {
                    loopTasks.submit(() -> {
                        internalTime.start();
                        processPartialConfNodeAsync(internalNode, internalNode.getConfSearchNode());
                        internalTime.stop();
                        return internalTime.getTimeS();
                    },(timeS) -> {
                        internalNode.markUpdated();
                        onFinishedInternal(timeS);
                    });
                }
            }
        }


        int maxNodes = 1000;
        if (leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int) Math.floor(leafTimeAverage / internalTimeAverage));
        System.out.println("Max nodes is now "+maxNodes);
        populateQueuesAsync(asyncQueue, maxNodes, 1);
        if(epsilonBound <= targetEpsilon)
            return;
        debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
        internalZ = state.internalZ.multiply(new BigDecimal(bias)); //MathTools.bigDivide(ZSums[0], new BigDecimal(Math.max(1,internalTimeAverage*internalNodes.size())), PartitionFunction.decimalPrecision);
        leafZ = state.leafZ; //MathTools.bigDivide(ZSums[1], new BigDecimal(Math.max(1,leafTimeAverage)), PartitionFunction.decimalPrecision);
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e, dscore %4.4e, denergy %4.4e", internalZ, leafZ, state.dScore, state.dEnergy));
        System.out.println("Number of ongoing background corrections: "+state.correctingLeaves);
        if((internalNodes.size() >= maxNodes && MathTools.isLessThan(internalZ, leafZ) && state.dScore > state.dEnergy)){
                //&& state.correctingLeaves < 10000) {
            numNodes = leafNodes.size();
            System.out.println("Processing "+numNodes+" leaf nodes...");
            int maxLeaves = 1;
            int numLeaves = 0;
            for(MARKStarNode leafNode: leafNodes) {
                if(correctedNodeAsync(leafNode, leafNode.getConfSearchNode())) {
                    System.out.println("Corrected node. Not minimizing it this round.");
                    continue;
                }
                if(numLeaves >= maxLeaves) {
                    break;
                }
                loopTasks.submit(() -> processFullConfNodeAsync(leafNode, leafNode.getConfSearchNode()),(leafResult) -> {
                    leafNode.markUpdated();
                    onFinishedLeaf(leafResult);
                });
                leafNode.markUpdated();
                numLeaves++;
                debugPrint("Processing Node: " + leafNode.getConfSearchNode().toString());
            }
            state.leafNodes.clear();
            state.leafZ = BigDecimal.ZERO;
            state.dScore *= 2;
            synchronized (this) {
                if (bias <1 )
                    bias = 1;
                else bias *=2;
                System.out.println("Biasing internal nodes by a factor of "+bias);
            }
        }
        else {
            numNodes = internalNodes.size();
            if(internalNodes.size() >1) {
                System.out.println("Processing " + numNodes + " internal nodes...");
                for (MARKStarNode internalNode : internalNodes) {
                    Stopwatch internalTime = new Stopwatch();
                    if (!MathTools.isGreaterThan(internalNode.getLowerBound(), BigDecimal.ZERO) &&
                            MathTools.isGreaterThan(
                                    MathTools.bigDivide(internalNode.getUpperBound(), rootNode.getUpperBound(),
                                            PartitionFunction.decimalPrecision),
                                    new BigDecimal(1 - targetEpsilon))
                            ) {
                        loopTasks.submit(() -> {
                            internalTime.start();
                            boundLowestBoundConfUnderNodeAsync(internalNode);
                            internalTime.stop();
                            return internalTime.getTimeS();
                        }, (timeS) -> {
                            internalNode.markUpdated();
                            onFinishedInternal(timeS);
                        });
                    } else {
                        loopTasks.submit(() -> {
                            internalTime.start();
                            processPartialConfNodeAsync(internalNode, internalNode.getConfSearchNode());
                            internalTime.stop();
                            return internalTime.getTimeS();
                        }, (timeS) -> {
                            internalNode.markUpdated();
                            onFinishedInternal(timeS);
                        });
                    }
                    internalNode.markUpdated();
                    if(!state.corrections.isEmpty()) {
                        RCTuple tuple = state.corrections.poll();
                        if(tuple != null)
                        correctionTasks.submit(() -> {
                            System.out.println("Submit correction for "+tuple.stringListing());
                            computeDifference(tuple, minimizingEcalc);
                            return null;
                        }, (ignored) -> {
                        });
                    }
                }
            }
            synchronized (this) {
                if (bias >1 )
                    bias = 1;
                else bias /=2;
                System.out.println("Biasing internal nodes by a factor of "+bias);
            }
            internalNodes.clear();
            state.internalZ = BigDecimal.ZERO;
            state.dEnergy *= 2;
                /*
            if(state.corrections.size() > 0)
            {
                            for (RCTuple tup : state.corrections) {
                                loopTasks.submit(() -> {
                                computeDifference(tup, minimizingEcalc);
                                            return null;
                                        },(ignored) -> {}

                                );
                            }
                state.corrections.clear();
            }
            */
        }
        if (epsilonBound <= targetEpsilon)
            return;
        loopCleanupAsync(loopWatch, numNodes);
    }

    private void onFinishedLeaf(LeafResult result) {
        synchronized(this) {
            ConfSearch.ScoredConf conf = result.conf;
            double energy = result.energy;
            Node node = result.node;
            double timeS = result.timeS;
            double newConfLower = result.newConfLower;
            double oldgscore = result.oldgscore;
            numConfsEnergied++;
            minList.set(conf.getAssignments().length-1,minList.get(conf.getAssignments().length-1)+1);
            recordReduction(conf.getScore(), energy);
            leafTimeSum += timeS;
            numLeavesScored++;
            leafTimeAverage = leafTimeSum/numLeavesScored;
            System.out.println("Processed leaf in "+timeS+" seconds.");
            profilePrint("Leaf time average is now "+leafTimeAverage);
            double delta = epsilonBound;
            updateBound();
            printMinimizationOutput(node, newConfLower, oldgscore);
            state.dEnergy = epsilonBound - delta;
            if(state.dEnergy == 0)
                state.dEnergy = state.dScore/10;
            System.out.println("dEnergy set to "+state.dEnergy);
        }
    }

    private void onFinishedInternal(Double timeS) {
        synchronized (this) {
            double curDelta = epsilonBound;
            updateBound();
            state.dScore = epsilonBound - curDelta;
            numInternalScored++;
            internalTimeSum += timeS;
            internalTimeAverage = internalTimeSum / numInternalScored;
        }
    }



    private boolean correctedNodeAsync(MARKStarNode curNode, Node node) {
        assert(curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.rigidScore);
            curNode.markUpdated();
            asyncQueue.add(curNode);
            return true;
        }
        return false;
    }

    private void populateQueuesAsync(BlockingQueue<MARKStarNode> queue,
                                     int maxNodes,
                                     int maxMinimizations){
        List<MARKStarNode> leftoverLeaves = new ArrayList<>();
        List<MARKStarNode> internalNodes = state.internalNodes;
        List<MARKStarNode> leafNodes = state.leafNodes;
        while(!queue.isEmpty() && (internalNodes.size() < maxNodes || leafNodes.size() < maxMinimizations)){
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore = node.getConfLowerBound() - node.gscore;
            double confCorrection = Math.min(correctgscore, node.rigidScore) + hscore;
            if(!node.isMinimized() && node.getConfLowerBound() - confCorrection > 1e-5) {
                recordCorrection(node.getConfLowerBound(), correctgscore - node.gscore);

                node.gscore = correctgscore;
                if (confCorrection > node.rigidScore) {
                    System.out.println("Overcorrected"+SimpleConfSpace.formatConfRCs(node.assignments)+": " + confCorrection + " > " + node.rigidScore);
                    node.gscore = node.rigidScore;
                    confCorrection = node.rigidScore + hscore;
                }
                node.setBoundsFromConfLowerAndUpper(confCorrection, node.getConfUpperBound());
                curNode.markUpdated();
                leftoverLeaves.add(curNode);
                continue;
            }

            BigDecimal diff = curNode.getUpperBound().subtract(curNode.getLowerBound());
            if (node.getLevel() < RCs.getNumPos() && internalNodes.size() < maxNodes) {
                if(internalNodes.size() < maxNodes) {
                    internalNodes.add(curNode);
                    state.internalZ = state.internalZ.add(diff);
                }
                else leftoverLeaves.add(curNode);
            }
            else if(shouldMinimize(node) && !correctedNode(leftoverLeaves, curNode, node)) {
                if(leafNodes.size() < maxMinimizations) {
                    leafNodes.add(curNode);
                    state.leafZ = state.leafZ.add(diff);
                }
                else
                    leftoverLeaves.add(curNode);
            }

        }
        queue.addAll(leftoverLeaves);
    }

    private void boundLowestBoundConfUnderNodeAsync(MARKStarNode startNode) {
        Comparator<MARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getConfSearchNode().getConfLowerBound());
        PriorityQueue<MARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
        drillQueue.add(startNode);

        List<MARKStarNode> newNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        while(!drillQueue.isEmpty()) {
            numNodes++;
            MARKStarNode curNode = drillQueue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);
            /*
            if(node.getLevel() + 1 == RCs.getNumPos())
                preminimizePartialAsync(startNode);
                */

            if (node.getLevel() < RCs.getNumPos()) {
                MARKStarNode nextNode = drillDown(newNodes, curNode, node);
                newNodes.remove(nextNode);
                drillQueue.add(nextNode);
            }
            else
                newNodes.add(curNode);

            //debugHeap(drillQueue, true);
            if(leafLoop.getTimeS() > 10) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",numNodes, overallLoop.getTime(2),rootNode.getLowerBound(),rootNode.getUpperBound()));
            }
        }
        asyncQueue.addAll(newNodes);

    }

    private void preminimizePartialAsync(MARKStarNode startNode) {
        if(MathTools.isGreaterThan(
                MathTools.bigDivide(startNode.getUpperBound(),rootNode.getUpperBound(),
                        PartitionFunction.decimalPrecision),
                new BigDecimal((1-targetEpsilon)/100000))) {
            RCTuple confTuple = startNode.toTuple();
            if(confTuple.size() +1 < RCs.getNumPos ())
                return;
            double threshhold = 0.5;
            double rigid = rigidEmat.getInternalEnergy(confTuple);
            double lower = minimizingEmat.getInternalEnergy(confTuple);
            if (rigid - lower > threshhold){
                computeTupleCorrection(minimizingEcalc, startNode.toTuple());
            }
        }
    }

    private void loopCleanupAsync(Stopwatch loopWatch, int numNodes) {
        loopWatch.stop();
        double loopTime = loopWatch.getTimeS();
        //profilePrint("Processed "+numNodes+" this loop.");
        loopWatch.reset();
        loopWatch.start();
        //processPreminimization(minimizingEcalc);
        //profilePrint("Preminimization time : "+loopWatch.getTime(2));
        double curEpsilon = epsilonBound;
        //rootNode.updateConfBounds(new ConfIndex(RCs.getNumPos()), RCs, gscorer, hscorer);
        loopWatch.stop();
        cleanupTime = loopWatch.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
        //System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
    }

    private void processPartialConfNodeAsync(MARKStarNode curNode, Node node) {
        // which pos to expand next?
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            ConfIndex confIndex = context.index;
            node.index(confIndex);
            int nextPos = order.getNextPos(confIndex, RCs);
            assert (!confIndex.isDefined(nextPos));
            assert (confIndex.isUndefined(nextPos));

            //preminimizePartialAsync(curNode);

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


                Stopwatch partialTime = new Stopwatch().start();
                node.index(context.index);
                Node child = node.assign(nextPos, nextRc);

                // score the child node differentially against the parent node
                if (child.getLevel() < RCs.getNumPos()) {
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    child.gscore = diff;
                    //Correct for incorrect gscore.
                    rigiddiff = rigiddiff - node.gscore + node.rigidScore;
                    child.rigidScore = rigiddiff;

                    double confLowerBound = child.gscore + hdiff;
                    double confUpperbound = rigiddiff + maxhdiff;
                    child.computeNumConformations(RCs);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    confRigid = confRigid - node.gscore + node.rigidScore;

                    child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    checkBounds(confCorrection, confRigid);
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    child.gscore = child.getConfLowerBound();
                    child.rigidScore = confRigid;
                    numConfsScored++;
                    progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                }
                partialTime.stop();
                loopPartialTime += partialTime.getTimeS();


                if (Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                // collect the possible children
                if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                    children.add(MARKStarNodeChild);
                }
                if (!child.isMinimized()) {
                    asyncQueue.add(MARKStarNodeChild);
                } else
                    MARKStarNodeChild.computeEpsilonErrorBounds();

            }
            curNode.markUpdated();
        }
    }

    private LeafResult processFullConfNodeAsync(MARKStarNode curNode, Node node) {
        double leafTime = 0;
        LeafResult result = new LeafResult();
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            node.index(context.index);

            Stopwatch leafTimer = new Stopwatch().start();
            ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound());
            ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
            Stopwatch correctionTimer = new Stopwatch().start();

            double energy = analysis.epmol.energy;
            double newConfUpper = energy;
            double newConfLower = energy;
            checkConfLowerBound(node, energy);
            if (energy > node.getConfUpperBound()) {
                System.err.println("Upper bounds got worse after minimization:" + energy
                        + " > " + (node.getConfUpperBound())+". Rejecting minimized energy.");
                System.err.println("Node info: "+node);

                newConfUpper = node.getConfUpperBound();
                newConfLower = node.getConfUpperBound();
            }
            curNode.setBoundsFromConfLowerAndUpper(newConfLower,newConfUpper);
            double oldgscore = node.gscore;
            node.gscore = newConfLower;
            String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound()) + "," + (node.getConfUpperBound()) + "]";
            debugPrint(out);
            leafTimer.stop();
            computeEnergyCorrection(analysis, conf);
            curNode.markUpdated();
            leafTime = leafTimer.getTimeS();
            result.conf = conf;
            result.node = node;
            result.newConfLower = newConfLower;
            result.oldgscore = oldgscore;
            result.energy = energy;
            result.timeS = leafTime;
        }
        // Dummy function. We're not doing anything here.
        progress.reportLeafNode(node.gscore, queue.size(), epsilonBound);
        if(!node.isMinimized())
            asyncQueue.add(curNode);

        return result;
    }

    private static class State {

        public Queue<RCTuple> corrections = new LinkedBlockingQueue<>();
        int correctingLeaves;
        BigDecimal numConfs;

		// upper bound (score axis) vars
		long numScoredConfs = 0;
		BigDecimal upperScoreWeightSum = BigDecimal.ZERO;
		BigDecimal minUpperScoreWeight = MathTools.BigPositiveInfinity;

		// lower bound (energy axis) vars
		long numEnergiedConfs = 0;
		BigDecimal lowerScoreWeightSum = BigDecimal.ZERO;
		BigDecimal energyWeightSum = BigDecimal.ZERO;
		BigDecimal minLowerScoreWeight = MathTools.BigPositiveInfinity;
		BigDecimal cumulativeZReduction = BigDecimal.ZERO;
		ArrayList<Integer> minList = new ArrayList<Integer>();
		List<MARKStarNode> internalNodes = new ArrayList<>();
		List<MARKStarNode> leafNodes = new ArrayList<>();
		BigDecimal internalZ = BigDecimal.ZERO;
		BigDecimal leafZ = BigDecimal.ZERO;

		// estimate of inital rates
		// (values here aren't super imporant since they get tuned during execution,
		// but make scores much faster than energies so we don't get a slow start)
		double scoreOps = 100.0;
		double energyOps = 1.0;

		double prevDelta = 1.0;
		double dEnergy = -1.0;
		double dScore = -1.0;

		State(BigInteger numConfs) {
			this.numConfs = new BigDecimal(numConfs);
		}

		double calcDelta() {
			BigDecimal upperBound = getUpperBound();
			if (MathTools.isZero(upperBound) || MathTools.isInf(upperBound)) {
				return 1.0;
			}
			return new BigMath(PartitionFunction.decimalPrecision)
				.set(upperBound)
				.sub(getLowerBound())
				.div(upperBound)
				.get()
				.doubleValue();
		}

		public BigDecimal getLowerBound() {
			return lowerScoreWeightSum;
		}

		public void printBoundStats() {
            System.out.println("Num confs: " + String.format("%12e",numConfs));
            System.out.println("Num Scored confs: " + String.format("%4d",numScoredConfs));
            String upperScoreString = minUpperScoreWeight.toString();
            String upperSumString = upperScoreWeightSum.toString();
            if(!MathTools.isInf(minUpperScoreWeight))
                upperScoreString = String.format("%12e",minUpperScoreWeight);
            if(!MathTools.isInf(upperScoreWeightSum))
                upperSumString = String.format("%12e",upperScoreWeightSum);
            System.out.println("Conf bound: " + upperScoreString);
            System.out.println("Scored weight bound:"+ upperSumString);
		}

		public BigDecimal getUpperBound() {
            return upperScoreWeightSum;
		}

		boolean epsilonReached(double targetEpsilon) {
			return calcDelta() <= targetEpsilon;
		}

		boolean isStable(BigDecimal stabilityThreshold) {
			return numEnergiedConfs <= 0 || stabilityThreshold == null || MathTools.isGreaterThanOrEqual(getUpperBound(), stabilityThreshold);
		}

		boolean hasLowEnergies() {
			return MathTools.isGreaterThan(minLowerScoreWeight,  BigDecimal.ZERO);
		}

		@Override
		public String toString() {
			return String.format("upper: count %d  sum %e  min %e     lower: count %d  score sum %e  energy sum %e",
				numScoredConfs, upperScoreWeightSum, minUpperScoreWeight,
				numEnergiedConfs, lowerScoreWeightSum, energyWeightSum
			);
		}
	}


    private boolean correctedNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        assert(curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.rigidScore);
            curNode.markUpdated();
            newNodes.add(curNode);
            return true;
        }
        return false;
    }

    private MARKStarNode drillDown(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        try(ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            ConfIndex confIndex = context.index;
            node.index(confIndex);
            // which pos to expand next?
            int nextPos = order.getNextPos(confIndex, RCs);
            assert (!confIndex.isDefined(nextPos));
            assert (confIndex.isUndefined(nextPos));

            // score child nodes with tasks (possibly in parallel)
            List<MARKStarNode> children = new ArrayList<>();
            double bestChildLower = Double.POSITIVE_INFINITY;
            MARKStarNode bestChild = null;
            for (int nextRc : RCs.get(nextPos)) {

                if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                    continue;
                }

                // if this child was pruned dynamically, then don't score it
                if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                    continue;
                }
                Stopwatch partialTime = new Stopwatch().start();
                Node child = node.assign(nextPos, nextRc);
                double confLowerBound = Double.POSITIVE_INFINITY;

                // score the child node differentially against the parent node
                if (child.getLevel() < RCs.getNumPos()) {
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double hdiff = context.hscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double maxhdiff = -context.negatedhscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    child.gscore = diff;
                    //Correct for incorrect gscore.
                    rigiddiff = rigiddiff - node.gscore + node.rigidScore;
                    child.rigidScore = rigiddiff;

                    confLowerBound = child.gscore + hdiff;
                    double confUpperbound = rigiddiff + maxhdiff;
                    child.computeNumConformations(RCs);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    confRigid = confRigid - node.gscore + node.rigidScore;

                    child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    checkBounds(confCorrection, confRigid);
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    child.gscore = child.getConfLowerBound();
                    confLowerBound = lowerbound;
                    child.rigidScore = confRigid;
                    numConfsScored++;
                    progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                }
                partialTime.stop();
                loopPartialTime += partialTime.getTimeS();


                if (Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                MARKStarNodeChild.markUpdated();
                if (confLowerBound < bestChildLower) {
                    bestChild = MARKStarNodeChild;
                }
                // collect the possible children
                if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                    children.add(MARKStarNodeChild);
                }
                newNodes.add(MARKStarNodeChild);

            }
            return bestChild;
        }
    }

    private void boundLowestBoundConfUnderNode(MARKStarNode startNode, List<MARKStarNode> generatedNodes) {
        Comparator<MARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getConfSearchNode().getConfLowerBound());
        PriorityQueue<MARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
        drillQueue.add(startNode);

        List<MARKStarNode> newNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        while(!drillQueue.isEmpty()) {
            numNodes++;
            MARKStarNode curNode = drillQueue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);

            if (node.getLevel() < RCs.getNumPos()) {
                MARKStarNode nextNode = drillDown(newNodes, curNode, node);
                newNodes.remove(nextNode);
                drillQueue.add(nextNode);
            }
            else
                newNodes.add(curNode);

            //debugHeap(drillQueue, true);
            if(leafLoop.getTimeS() > 10) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",numNodes, overallLoop.getTime(2),rootNode.getLowerBound(),rootNode.getUpperBound()));
            }
        }
        generatedNodes.addAll(newNodes);

    }


    private static class CorrectionData {
        ConfAnalyzer.ConfAnalysis analysis;
        ConfSearch.ScoredConf conf;
    }



    private void printMinimizationOutput(Node node, double newConfLower, double oldgscore) {
        if (printMinimizedConfs) {
            System.out.println("[" + SimpleConfSpace.formatConfRCs(node.assignments) + "]"
                    + String.format("conf:%4d, score:%12.6f, lower:%12.6f, corrected:%12.6f energy:%12.6f"
                            +", bounds:[%12e, %12e], delta:%12.6f, time:%10s",
                    numConfsEnergied, oldgscore, minimizingEmat.confE(node.assignments),
                    correctionMatrix.confE(node.assignments), newConfLower,
                    rootNode.getConfSearchNode().getSubtreeLowerBound(),rootNode.getConfSearchNode().getSubtreeUpperBound(),
                    epsilonBound, stopwatch.getTime(2)));

        }
    }

    private void checkConfLowerBound(Node node, double energy) {
        if(energy < node.getConfLowerBound()) {
            System.err.println("Bounds are incorrect:" + (node.getConfLowerBound()) + " > "
                    + energy);
            if (energy < 10)
                System.err.println("The bounds are probably wrong.");
            //System.exit(-1);
        }
    }


    private void checkBounds(double lower, double upper)
    {
        if (upper < lower && upper - lower > 1e-5 && upper < 10)
            debugPrint("Bounds incorrect.");
    }

    private void computeEnergyCorrection(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf){
        if(conf.getAssignments().length < 3)
            return;
        // TODO: Replace the sortedPairwiseTerms with an ArrayList<TupE>.
        //System.out.println("Analysis:"+analysis);
        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(minimizingEmat);
        Stopwatch correctionTime = new Stopwatch().start();
        //System.out.println("Energy Analysis: "+energyAnalysis);
        //System.out.println("Score Analysis: "+scoreAnalysis);
        EnergyMatrix diff = energyAnalysis.diff(scoreAnalysis);
        //System.out.println("Difference Analysis " + diff);
        List<TupE> sortedPairwiseTerms2 = new ArrayList<>();
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
                        TupE tupe = new TupE(new RCTuple(pos, rc, pos2, rc2), sum);
                        sortedPairwiseTerms2.add(tupe);
                    }
                }
            }
        }
        Collections.sort(sortedPairwiseTerms2);

        double threshhold = 0.1;
        double minDifference = 0.9;
        double maxDiff = sortedPairwiseTerms2.get(0).E;
        for(int i = 0; i < sortedPairwiseTerms2.size(); i++)
        {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if(pairDiff < minDifference || maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            int localMinimizations = 0;
            for(int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double tupleBounds = rigidEmat.getInternalEnergy(tuple) - minimizingEmat.getInternalEnergy(tuple);
                if(tupleBounds < minDifference)
                    continue;
                computeDifference(tuple, minimizingEcalc);
                minList.set(tuple.size()-1,minList.get(tuple.size()-1)+1);
                localMinimizations++;
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        correctionTime.stop();
        //ecalc.tasks.waitForFinish();
    }




    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if(correctedTuples.contains(tuple.stringListing()))
            return;
        correctedTuples.add(tuple.stringListing());
        if(correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        synchronized (this) {
            state.correctingLeaves++;
        }
        double tripleEnergy = minimizingEcalc.calcEnergy(tuple).energy;

        double lowerbound = minimizingEmat.getInternalEnergy(tuple);
        if (tripleEnergy - lowerbound > 0) {
            double correction = tripleEnergy - lowerbound;
            correctionMatrix.setHigherOrder(tuple, correction);
        }
        else
            System.err.println("Negative correction for "+tuple.stringListing());
        synchronized (this) {
            state.correctingLeaves--;
        }
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for(int pos: positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    private void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap) {
        if(correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        System.out.println("Computing correction for " + overlap.stringListing() + " penalty of " + (partiallyMinimizedLower - pairwiseLower));
        progress.reportPartialMinimization(1, epsilonBound);
        if(partiallyMinimizedLower > pairwiseLower)
        synchronized (correctionMatrix) {
            correctionMatrix.setHigherOrder(overlap, partiallyMinimizedLower - pairwiseLower);
        }
        progress.reportPartialMinimization(1, epsilonBound);
    }

    private List<MARKStarNode> getTopConfs(int numConfs) {
        List<MARKStarNode> topConfs = new ArrayList<>();
        while (topConfs.size() < numConfs&& !queue.isEmpty()) {
            MARKStarNode nextLowestConf = queue.poll();
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

    private void updateBound() {
        double curEpsilon = epsilonBound;
        Stopwatch time = new Stopwatch().start();
        epsilonBound = rootNode.computeEpsilonErrorBounds();
        time.stop();
        //profilePrint("Bound update time: "+time.getTimeS());
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

    private class LeafResult {
        public ConfSearch.ScoredConf conf;
        public double oldgscore;
        public double energy;
        public Node node;
        public double newConfLower;
        public double timeS;
    }
}
