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
import edu.duke.cs.osprey.astar.conf.order.DynamicHMeanAStarOrder;
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
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;

public class MARKStarBound implements PartitionFunction {

    protected double targetEpsilon = 1;
    public boolean debug = false;
    public boolean profileOutput = false;
    private Status status = null;
    private Values values = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;
    // max confs minimized, -1 means infinite.
    private int maxNumConfs = -1;

    protected int maxMinimizations = 1;


    // the number of full conformations scored OR energied
    private int numConfsScored = 0;

    protected int numInternalNodesProcessed = 0;

    private boolean printMinimizedConfs;
    private MARKStarProgress progress;
    public String stateName = String.format("%4f",Math.random());
    private int numPartialMinimizations;
    private ArrayList<Integer> minList;
    protected double internalTimeAverage;
    protected double leafTimeAverage;
    private double cleanupTime;
    private boolean nonZeroLower;
    protected static TaskExecutor loopTasks;

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
    public void setStabilityThreshold(BigDecimal threshold) {
        stabilityThreshold = threshold;
    }

    public void setMaxNumConfs(int maxNumConfs) {
        this.maxNumConfs = maxNumConfs;
    }

    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        status = Status.Estimating;
        values = new Values();
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        targetEpsilon = epsilon;
        status = Status.Estimating;
        values = new Values();
        this.stabilityThreshold = stabilityThreshold;
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

    private int workDone() {
        return numInternalNodesProcessed + numConfsEnergied + numConfsScored + numPartialMinimizations ;
    }

    @Override
    public void compute(int maxNumConfs) {
        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;

        int previousConfCount = workDone();

        if(!nonZeroLower) {
            runUntilNonZero();
            updateBound();
        }
        while (epsilonBound > targetEpsilon &&
                workDone()-previousConfCount < maxNumConfs
                && isStable(stabilityThreshold)) {
            debugPrint("Tightening from epsilon of "+epsilonBound);
            if(debug)
                debugHeap(queue);
            tightenBoundInPhases();
            debugPrint("Errorbound is now "+epsilonBound);
            if(lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("Error. Bounds got looser.");
                //System.exit(-1);
            }
            lastEps = epsilonBound;
        }
        if(!isStable(stabilityThreshold))
            status = Status.Unstable;
        loopTasks.waitForFinish();
        minimizingEcalc.tasks.waitForFinish();
        BigDecimal averageReduction = BigDecimal.ZERO;
        int totalMinimizations = numConfsEnergied + numPartialMinimizations;
        if(totalMinimizations> 0)
            averageReduction = cumulativeZCorrection
                .divide(new BigDecimal(totalMinimizations), new MathContext(BigDecimal.ROUND_HALF_UP));
        debugPrint(String.format("Average Z reduction per minimization: %12.6e",averageReduction));
        values.pstar = rootNode.getUpperBound();
        values.qstar = rootNode.getLowerBound();
        values.qprime= rootNode.getUpperBound();
        if(epsilonBound < targetEpsilon) {
            status = Status.Estimated;
            if(values.qstar.compareTo(BigDecimal.ZERO) == 0) {
                status = Status.Unstable;
            }
            //rootNode.printTree(stateName, minimizingEcalc.confSpace);
        }
    }

    protected void debugPrint(String s) {
        if(debug)
            System.out.println(s);
    }

    protected void profilePrint(String s) {
        if(profileOutput)
            System.out.println(s);
    }

    public void compute() {
        compute(Integer.MAX_VALUE);
    }

    @Override
    public Result makeResult() {
        // Calculate the upper bound z reductions from conf lower bounds, since we don't explicitly record these
        lowerReduction_ConfUpperBound = rootNode.getLowerBound().subtract(startLowerBound).subtract(lowerReduction_FullMin);
        // Calculate the lower bound z reductions from conf upper bounds, since we don't explicitly record these
        upperReduction_ConfLowerBound = startUpperBound.subtract(rootNode.getUpperBound()).subtract(upperReduction_FullMin).subtract(upperReduction_PartialMin);

        PartitionFunction.Result result = new PartitionFunction.Result(getStatus(), getValues(), getNumConfsEvaluated());
        /*
        result.setWorkInfo(numPartialMinimizations, numConfsScored,minList);
        result.setZInfo(lowerReduction_FullMin, lowerReduction_ConfUpperBound, upperReduction_FullMin, upperReduction_PartialMin, upperReduction_ConfLowerBound);
        result.setOrigBounds(startUpperBound, startLowerBound);
        result.setTimeInfo(stopwatch.getTimeNs());
        result.setMiscInfo(new BigDecimal(rootNode.getNumConfs()));
        */
        return result;
    }


    /**
     * TODO: 1. Make MARKStarBounds use and update a queue.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */
    // We keep track of the root node for computing our K* bounds
    protected MARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    protected final Queue<MARKStarNode> queue;
    protected double epsilonBound = Double.POSITIVE_INFINITY;
    private ConfIndex confIndex;
    public final AStarOrder order;
    // TODO: Implement new AStarPruner for MARK*?
    public final AStarPruner pruner;
    protected RCs RCs;
    protected Parallelism parallelism;
    private ObjectPool<ScoreContext> contexts;
    private MARKStarNode.ScorerFactory gscorerFactory;
    private MARKStarNode.ScorerFactory hscorerFactory;

    public boolean reduceMinimizations = true;
    private ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix rigidEmat;
    UpdatingEnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc;
    private Stopwatch stopwatch = new Stopwatch().start();
    // Variables for reporting pfunc reductions more accurately
    BigDecimal startUpperBound = null; //can't start with infinity
    BigDecimal startLowerBound = BigDecimal.ZERO;
    BigDecimal lowerReduction_FullMin = BigDecimal.ZERO; //Pfunc lower bound improvement from full minimization
    BigDecimal lowerReduction_ConfUpperBound = BigDecimal.ZERO; //Pfunc lower bound improvement from conf upper bounds
    BigDecimal upperReduction_FullMin = BigDecimal.ZERO; //Pfunc upper bound improvement from full minimization
    BigDecimal upperReduction_PartialMin = BigDecimal.ZERO; //Pfunc upper bound improvement from partial minimization corrections
    BigDecimal upperReduction_ConfLowerBound = BigDecimal.ZERO; //Pfunc upper bound improvement from conf lower bounds

    BigDecimal cumulativeZCorrection = BigDecimal.ZERO;//Pfunc upper bound improvement from partial minimization corrections
    BigDecimal ZReductionFromMin = BigDecimal.ZERO;//Pfunc lower bound improvement from full minimization
    BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private boolean computedCorrections = false;
    private long loopPartialTime = 0;
    private Set<String> correctedTuples = Collections.synchronizedSet(new HashSet<>());
    private BigDecimal stabilityThreshold;
    private double leafTimeSum = 0;
    private double internalTimeSum = 0;
    private int numLeavesScored = 0;
    private int numInternalScored = 0;

    public static MARKStarBound makeFromConfSpaceInfo(BBKStar.ConfSpaceInfo info, RCs rcs) {
        throw new UnsupportedOperationException("MARK* is not yet integrated into BBK*. Coming soon!");
        /*
        ConfEnergyCalculator minimizingConfEcalc = info.confEcalcMinimized;
        return new MARKStarBound(info.confSpace, info.ematRigid, info.ematMinimized, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
        */
    }

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.queue = new PriorityQueue<>();
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        MPLPPairwiseHScorer scorer = new MPLPPairwiseHScorer(updater, minimizingEmat, 1, 0.0001);
        hscorerFactory = (emats) -> new MPLPPairwiseHScorer(updater, emats, 1, 0.0001);//TraditionalPairwiseHScorer(emats, rcs);//

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs,
                gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat),
                gscorerFactory.make(rigidEmat),
                new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);
                //hscorerFactory.make(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs), true);
        confIndex = new ConfIndex(rcs.getNumPos());
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.RCs = rcs;
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat),hscorerFactory.make(minimizingEmat));
        this.pruner = null;

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.gscorer = gscorerFactory.make(minimizingEmat);
            context.hscorer = hscorerFactory.make(minimizingEmat);
            context.rigidscorer = gscorerFactory.make(rigidEmat);
            /** These scoreres should match the scorers in the MARKStarNode root - they perform the same calculations**/
            context.negatedhscorer = new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        progress = new MARKStarProgress(RCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        setParallelism(parallelism);
        updateBound();
        // Recording pfunc starting bounds
        this.startLowerBound = rootNode.getLowerBound();
        this.startUpperBound = rootNode.getUpperBound();
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(),0));
    }

    protected static class ScoreContext {
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
        //loopTasks = minimizingEcalc.tasks;
        if(loopTasks == null)
            loopTasks = parallelism.makeTaskExecutor(1000);
        contexts.allocate(parallelism.getParallelism());
    }

    private void debugEpsilon(double curEpsilon) {
        if(debug && curEpsilon < epsilonBound) {
            System.err.println("Epsilon just got bigger.");
        }
    }

    protected boolean shouldMinimize(Node node) {
        return node.getLevel() == RCs.getNumPos() && !node.isMinimized();
    }

    protected void recordCorrection(double lowerBound, double correction) {
        BigDecimal upper = bc.calc(lowerBound);
        BigDecimal corrected = bc.calc(lowerBound + correction);
        cumulativeZCorrection = cumulativeZCorrection.add(upper.subtract(corrected));
        upperReduction_PartialMin = upperReduction_PartialMin.add(upper.subtract(corrected));
    }
    private void recordReduction(double lowerBound, double upperBound, double energy) {
        BigDecimal lowerBoundWeight = bc.calc(lowerBound);
        BigDecimal upperBoundWeight = bc.calc(upperBound);
        BigDecimal energyWeight = bc.calc(energy);
        ZReductionFromMin = ZReductionFromMin.add(lowerBoundWeight.subtract(upperBoundWeight));
        upperReduction_FullMin = upperReduction_FullMin.add(lowerBoundWeight.subtract(energyWeight));
        lowerReduction_FullMin = lowerReduction_FullMin.add(energyWeight.subtract(upperBoundWeight));

    }

    private void debugBreakOnConf(int[] conf) {
        int[] confOfInterest = new int[]{4,5,8,18};
        if(conf.length != confOfInterest.length)
            return;
        boolean match = true;
        for(int i = 0; i < confOfInterest.length; i++) {
            if(conf[i] != confOfInterest[i]) {
                match = false;
                break;
            }
        }
        if(match)
            System.out.println("Matched "+SimpleConfSpace.formatConfRCs(conf));
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


        newNodes.clear();
        System.out.println("Found a leaf!");
        nonZeroLower = true;
    }

    protected void tightenBoundInPhases() {
        System.out.println(String.format("Current overall error bound: %12.10f, spread of [%12.6e, %12.6e]",epsilonBound, rootNode.getLowerBound(), rootNode.getUpperBound()));
        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        int numNodes = 0;
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        Stopwatch internalTime = new Stopwatch();
        Stopwatch leafTime = new Stopwatch();
        double leafTimeSum = 0;
        double internalTimeSum = 0;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ,leafZ};
        populateQueues(queue, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        updateBound();
        debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
        internalZ = ZSums[0];// MathTools.bigDivide(ZSums[0], new BigDecimal(Math.max(1,internalTimeAverage*internalNodes.size())), PartitionFunction.decimalPrecision);
        leafZ = ZSums[1]; //MathTools.bigDivide(ZSums[1], new BigDecimal(Math.max(1,leafTimeAverage)), PartitionFunction.decimalPrecision);
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        if(MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            System.out.println("Processing "+numNodes+" leaf nodes...");
            leafTime.reset();
            leafTime.start();
            for(MARKStarNode leafNode: leafNodes) {
                processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
                debugPrint("Processing Node: " + leafNode.getConfSearchNode().toString());
            }
            loopTasks.waitForFinish();
            leafTime.stop();
            leafTimeAverage = leafTime.getTimeS();
            System.out.println("Processed "+numNodes+" leaves in "+leafTimeAverage+" seconds.");
            queue.addAll(internalNodes);
            if(maxMinimizations < parallelism.numThreads)
                maxMinimizations++;
        }
        else {
            numNodes = internalNodes.size();
            System.out.println("Processing "+numNodes+" internal nodes...");
            internalTime.reset();
            internalTime.start();
            for (MARKStarNode internalNode : internalNodes) {
                if(!MathTools.isGreaterThan(internalNode.getLowerBound(),BigDecimal.ONE) &&
                    MathTools.isGreaterThan(
                            MathTools.bigDivide(internalNode.getUpperBound(),rootNode.getUpperBound(),
                                    PartitionFunction.decimalPrecision),
                            new BigDecimal(1-targetEpsilon))
                ) {
                    loopTasks.submit(() -> {
                        boundLowestBoundConfUnderNode(internalNode, newNodes);
                        return null;
                    }, (ignored) -> {
                    });
                }
                else {
                    processPartialConfNode(newNodes, internalNode, internalNode.getConfSearchNode());
                }
                internalNode.markUpdated();
            }
            loopTasks.waitForFinish();
            internalTime.stop();
            internalTimeSum=internalTime.getTimeS();
            internalTimeAverage = internalTimeSum/Math.max(1,internalNodes.size());
            debugPrint("Internal node time :"+internalTimeSum+", average "+internalTimeAverage);
            queue.addAll(leafNodes);
            numInternalNodesProcessed+=internalNodes.size();
        }
        if (epsilonBound <= targetEpsilon)
            return;
        loopCleanup(newNodes, loopWatch, numNodes);
    }

    protected void debugHeap(Queue<MARKStarNode> queue) {
        int maxNodes = 10;
        System.out.println("Node heap:");
        List<MARKStarNode> nodes = new ArrayList<>();
        while(!queue.isEmpty() && nodes.size() < 10)
        {
            MARKStarNode node = queue.poll();
            System.out.println(node.getConfSearchNode());
            nodes.add(node);
        }
        queue.addAll(nodes);
    }


    boolean isStable(BigDecimal stabilityThreshold) {
        return numConfsEnergied <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(rootNode.getUpperBound(), stabilityThreshold);
    }


    protected void populateQueues(Queue<MARKStarNode> queue, List<MARKStarNode> internalNodes, List<MARKStarNode> leafNodes, BigDecimal internalZ,
                                BigDecimal leafZ, BigDecimal[] ZSums) {
        List<MARKStarNode> leftoverLeaves = new ArrayList<>();
        int maxNodes = 1000;
        if(leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int)Math.floor(0.1*leafTimeAverage/internalTimeAverage));
        while(!queue.isEmpty() && internalNodes.size() < maxNodes){
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore = node.getConfLowerBound() - node.gscore;
            double confCorrection = Math.min(correctgscore, node.rigidScore) + hscore;
            if(!node.isMinimized() && node.getConfLowerBound() < confCorrection
                    && node.getConfLowerBound() - confCorrection > 1e-5) {
                if(confCorrection < node.getConfLowerBound()) {
                    System.out.println("huh!?");
                }
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
            if (node.getLevel() < RCs.getNumPos()) {
                internalNodes.add(curNode);
                internalZ = internalZ.add(diff);
            }
            else if(shouldMinimize(node) && !correctedNode(leftoverLeaves, curNode, node)) {
                if(leafNodes.size() < maxMinimizations) {
                    leafNodes.add(curNode);
                    leafZ = leafZ.add(diff);
                }
                else
                    leftoverLeaves.add(curNode);
            }

        }
        ZSums[0] = internalZ;
        ZSums[1] = leafZ;
        queue.addAll(leftoverLeaves);
    }

    protected void loopCleanup(List<MARKStarNode> newNodes, Stopwatch loopWatch, int numNodes) {
        for(MARKStarNode node: newNodes) {
            if(node != null)
                queue.add(node);
        }
        loopWatch.stop();
        double loopTime = loopWatch.getTimeS();
        profilePrint("Processed "+numNodes+" this loop, spawning "+newNodes.size()+" in "+loopTime+", "+stopwatch.getTime()+" so far");
        loopWatch.reset();
        loopWatch.start();
        processPreminimization(minimizingEcalc);
        profilePrint("Preminimization time : "+loopWatch.getTime(2));
        double curEpsilon = epsilonBound;
        //rootNode.updateConfBounds(new ConfIndex(RCs.getNumPos()), RCs, gscorer, hscorer);
        updateBound();
        loopWatch.stop();
        cleanupTime = loopWatch.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
        System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
    }

    protected boolean correctedNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        assert(curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if((node.getLevel() == RCs.getNumPos() && node.getConfLowerBound()< confCorrection)
                || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            node.setBoundsFromConfLowerAndUpper(node.getConfLowerBound() - oldg + confCorrection, node.getConfUpperBound());
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
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                    progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    confRigid=confRigid-node.gscore+node.rigidScore;

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

    protected void boundLowestBoundConfUnderNode(MARKStarNode startNode, List<MARKStarNode> generatedNodes) {
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
            else {
                newNodes.add(curNode);
            }

            //debugHeap(drillQueue, true);
            if(leafLoop.getTimeS() > 1) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",numNodes, overallLoop.getTime(2),rootNode.getLowerBound(),rootNode.getUpperBound()));
            }
        }
        generatedNodes.addAll(newNodes);

    }

    protected void processPartialConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        // which pos to expand next?
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

            loopTasks.submit(() -> {

                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    Stopwatch partialTime = new Stopwatch().start();
                    ScoreContext context = checkout.get();
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
                        rigiddiff=rigiddiff-node.gscore+node.rigidScore;
                        child.rigidScore = rigiddiff;

                        double confLowerBound = child.gscore + hdiff;
                        double confUpperbound = rigiddiff + maxhdiff;
                        child.computeNumConformations(RCs);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if(diff < confCorrection) {
                            recordCorrection(confLowerBound, confCorrection - diff);
                            confLowerBound = confCorrection + hdiff;
                        }
                        child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                        progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        confRigid=confRigid-node.gscore+node.rigidScore;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = minimizingEmat.confE(child.assignments);

                        if(lowerbound < confCorrection) {
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection,confRigid);
                        child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                        child.gscore = confCorrection;
                        child.rigidScore = confRigid;
                        numConfsScored++;
                        progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                    }
                    partialTime.stop();
                    loopPartialTime+=partialTime.getTimeS();


                    return child;
                }

            }, (Node child) -> {
                if(Double.isNaN(child.rigidScore))
                    System.out.println("Huh!?");
                MARKStarNode MARKStarNodeChild = curNode.makeChild(child);
                    // collect the possible children
                    if (MARKStarNodeChild.getConfSearchNode().getConfLowerBound() < 0) {
                        children.add(MARKStarNodeChild);
                    }
                    if (!child.isMinimized()) {
                        newNodes.add(MARKStarNodeChild);
                    }
                    else
                        MARKStarNodeChild.computeEpsilonErrorBounds();

                curNode.markUpdated();
            });
        }
    }


    protected void processFullConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        double confCorrection = correctionMatrix.confE(node.assignments);
        if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.getConfUpperBound());
            curNode.markUpdated();
            newNodes.add(curNode);
            return;
        }
        loopTasks.submit(() -> {
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                node.index(context.index);

                ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound());
                ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                Stopwatch correctionTimer = new Stopwatch().start();
                computeEnergyCorrection(analysis, conf, context.ecalc);

                double energy = analysis.epmol.energy;
                double newConfUpper = energy;
                double newConfLower = energy;
                // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
                double oldConfUpper = node.getConfUpperBound();
                double oldConfLower = node.getConfLowerBound();
                checkConfLowerBound(node, energy);
                if (newConfUpper > oldConfUpper) {
                    System.err.println("Upper bounds got worse after minimization:" + newConfUpper
                            + " > " + (oldConfUpper)+". Rejecting minimized energy.");
                    System.err.println("Node info: "+node);

                    newConfUpper = oldConfUpper;
                    newConfLower = oldConfUpper;
                }
                curNode.setBoundsFromConfLowerAndUpper(newConfLower,newConfUpper);
                double oldgscore = node.gscore;
                node.gscore = newConfLower;
                String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound()) + "," + (node.getConfUpperBound()) + "]";
                debugPrint(out);
                curNode.markUpdated();
                synchronized(this) {
                    numConfsEnergied++;
                    minList.set(conf.getAssignments().length-1,minList.get(conf.getAssignments().length-1)+1);
                    recordReduction(oldConfLower, oldConfUpper, energy);
                    printMinimizationOutput(node, newConfLower, oldgscore);
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

    private void computeEnergyCorrection(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                                  ConfEnergyCalculator ecalc) {
        if(conf.getAssignments().length < 3)
            return;
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
        double triplethreshhold = 0.3;
        double maxDiff = sortedPairwiseTerms2.get(0).E;
        for(int i = 0; i < sortedPairwiseTerms2.size(); i++)
        {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if(pairDiff < minDifference &&  maxDiff - pairDiff > threshhold)
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
                if(tupleBounds < triplethreshhold)
                    continue;
                minList.set(tuple.size()-1,minList.get(tuple.size()-1)+1);
                computeDifference(tuple, minimizingEcalc);
                localMinimizations++;
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        correctionTime.stop();
        ecalc.tasks.waitForFinish();
    }




    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if(correctedTuples.contains(tuple.stringListing()))
            return;
        correctedTuples.add(tuple.stringListing());
        if(correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
            double tripleEnergy = minimizedTuple.energy;

            double lowerbound = minimizingEmat.getInternalEnergy(tuple);
            if (tripleEnergy - lowerbound > 0) {
                double correction = tripleEnergy - lowerbound;
                correctionMatrix.setHigherOrder(tuple, correction);
            }
            else
                System.err.println("Negative correction for "+tuple.stringListing());
        });
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for(int pos: positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    private void processPreminimization(ConfEnergyCalculator ecalc) {
        int maxMinimizations = 1;//parallelism.numThreads;
        List<MARKStarNode> topConfs = getTopConfs(maxMinimizations);
        // Need at least two confs to do any partial preminimization
        if (topConfs.size() < 2) {
            queue.addAll(topConfs);
            return;
        }
        RCTuple lowestBoundTuple = topConfs.get(0).toTuple();
        RCTuple overlap = findLargestOverlap(lowestBoundTuple, topConfs, 3);
        //Only continue if we have something to minimize
        for (MARKStarNode conf : topConfs) {
            RCTuple confTuple = conf.toTuple();
            if(minimizingEmat.getInternalEnergy(confTuple) == rigidEmat.getInternalEnergy(confTuple))
                continue;
            numPartialMinimizations++;
            minList.set(confTuple.size()-1,minList.get(confTuple.size()-1)+1);
            if (confTuple.size() > 2 && confTuple.size() < RCs.getNumPos ()){
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(minimizingEcalc, conf.toTuple());
                    return null;
                }, (econf) -> {
                });
            }
        }
        //minimizingEcalc.tasks.waitForFinish();
        ConfIndex index = new ConfIndex(RCs.getNumPos());
        if(overlap.size() > 3 && !correctionMatrix.hasHigherOrderTermFor(overlap)
                && minimizingEmat.getInternalEnergy(overlap) != rigidEmat.getInternalEnergy(overlap)) {
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(ecalc, overlap);
                    return null;
                }, (econf) -> {
                });
        }
        queue.addAll(topConfs);
    }

    private void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap) {
        if(correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
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

    protected void updateBound() {
        double curEpsilon = epsilonBound;
        Stopwatch time = new Stopwatch().start();
        epsilonBound = rootNode.computeEpsilonErrorBounds();
        time.stop();
        //System.out.println("Bound update time: "+time.getTime(2));
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
