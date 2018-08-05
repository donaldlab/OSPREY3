package edu.duke.cs.osprey.markstar.framework;

import EDU.oswego.cs.dl.util.concurrent.FJTask;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.UpperLowerAStarOrder;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.TraditionalPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.EdgeUpdater;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.MPLPUpdater;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
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
import edu.duke.cs.osprey.tools.TimeFormatter;
import javafx.util.Pair;
import org.ojalgo.matrix.transformation.Rotation;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;

public class MARKStarBound implements PartitionFunction {

    private double targetEpsilon = 1;
    public boolean debug = false;
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
    private boolean hitLeaf;

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

    @Override
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
    public int getNumConfsScored() {
        return numConfsScored;
    }

    @Override
    public void compute(int maxNumConfs) {
        debugPrint("Num conformations: "+rootNode.getConfSearchNode().getNumConformations());
        double lastEps = 1;

        int previousConfCount = numConfsEnergied + numConfsScored + numPartialMinimizations;
        if(!hitLeaf) {
            runUntilLeaf();
            updateBound();
        }
        while (epsilonBound > targetEpsilon &&
                (maxNumConfs < 0 || numConfsEnergied + numConfsScored + numPartialMinimizations - previousConfCount < maxNumConfs)) {
            debugPrint("Tightening from epsilon of "+epsilonBound);
            tightenBoundInPhases();
            debugPrint("Errorbound is now "+epsilonBound);
            if(lastEps < epsilonBound && epsilonBound - lastEps > 0.01) {
                System.err.println("Error. Bounds got looser.");
                //System.exit(-1);
            }
            lastEps = epsilonBound;
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

    public void compute() {
        compute(-1);
    }

    @Override
    public Result makeResult() {
        Result result = new Result(getStatus(), getValues(), getNumConfsEvaluated(),numPartialMinimizations, numConfsScored, rootNode.getNumConfs(), Long.toString(stopwatch.getTimeNs()), minList, cumulativeZCorrection);
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
    private ObjectPool<ScoreContext> contexts;
    private MARKStarNode.ScorerFactory gscorerFactory;
    private MARKStarNode.ScorerFactory hscorerFactory;

    public static final int MAX_CONFSPACE_FRACTION = 1000000;
    public static final double MINIMIZATION_FACTOR = 0.1;
    public boolean reduceMinimizations = false;
    private ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix rigidEmat;
    UpdatingEnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc = null;
    private Stopwatch stopwatch = new Stopwatch().start();
    BigDecimal cumulativeZCorrection = BigDecimal.ZERO;
    BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private boolean computedCorrections = false;
    private long loopPartialTime = 0;

    public static MARKStarBound makeFromConfSpaceInfo(BBKStar.ConfSpaceInfo info, RCs rcs) {
        ConfEnergyCalculator minimizingConfEcalc = info.confEcalcMinimized;
        return new MARKStarBound(info.confSpace, info.ematRigid, info.ematMinimized, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
    }

    public MARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                         ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.queue = new PriorityQueue<>();
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);

        MPLPUpdater updater = new EdgeUpdater();
        hscorerFactory = (emats) -> new TraditionalPairwiseHScorer(emats, rcs);//MPLPPairwiseHScorer(updater, emats, 50, 0.03);

        rootNode = MARKStarNode.makeRoot(confSpace, rigidEmat, minimizingEmat, rcs, gscorerFactory, hscorerFactory, true);
        queue.add(rootNode);
        updateBound();
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
        leafTasks = parallelism.makeTaskExecutor(1000);
        internalTasks = parallelism.makeTaskExecutor(1000);
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

    // We want to process internal nodes without worrying about the bound too much until we have
    // a nonzero lower bound. We have to have a nonzero lower bound, so we have to have at least
    // one node with a negative conf upper bound.
    private void runUntilLeaf() {
        System.out.println("Running until leaf is found...");
        double bestConfUpper = Double.POSITIVE_INFINITY;

        List<MARKStarNode> newNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        int numNodes = 0;
        while(!queue.isEmpty() && leafNodes.isEmpty()){
            numNodes++;
            assert(newNodes.size() < 1);
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);

            bestConfUpper = Math.min(node.getConfUpperBound(), bestConfUpper);
            if (node.getLevel() < RCs.getNumPos()) {
                processPartialConfNode(newNodes, curNode, node);
            }
            else
                leafNodes.add(curNode);
            internalTasks.waitForFinish();
            queue.addAll(newNodes);
            newNodes.clear();
            if(numNodes % 1000 == 0)
                System.out.println("Processed "+numNodes+" nodes so far...");
        }
        queue.addAll(leafNodes);
        System.out.println("Found a leaf!");
        hitLeaf = true;
    }
    private void tightenBoundInPhases() {
        System.out.println(String.format("Current overall error bound: %12.6f",epsilonBound));
        List<MARKStarNode> internalNodes = new ArrayList<>();
        List<MARKStarNode> leafNodes = new ArrayList<>();
        List<MARKStarNode> newNodes = new ArrayList<>();
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
        populateQueues(internalNodes, leafNodes, internalZ, leafZ, ZSums);
        internalZ = MathTools.bigDivide(ZSums[0], new BigDecimal(Math.max(1,internalTimeAverage*internalNodes.size())), PartitionFunction.decimalPrecision);
        leafZ = MathTools.bigDivide(ZSums[1], new BigDecimal(Math.max(1,leafTimeAverage)), PartitionFunction.decimalPrecision);
        if(MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            leafTime.reset();
            leafTime.start();
            for(MARKStarNode leafNode: leafNodes) {
                processFullConfNode(newNodes, leafNode, leafNode.getConfSearchNode());
                debugPrint("Processing Node: " + leafNode.getConfSearchNode().toString());
            }
            leafTime.stop();
            leafTasks.waitForFinish();
            leafTimeSum = leafTime.getTimeS();
            leafTimeAverage = leafTimeSum/leafNodes.size();
            queue.addAll(internalNodes);
        }
        else {
            numNodes = internalNodes.size();
            internalTime.reset();
            internalTime.start();
            for (MARKStarNode internalNode : internalNodes) {
                processPartialConfNode(newNodes, internalNode, internalNode.getConfSearchNode());
                //debugPrint("Processing Node: " + internalNode.getConfSearchNode().toString());
            }
            internalTasks.waitForFinish();
            internalTime.stop();
            internalTimeSum=internalTime.getTimeS();
            internalTimeAverage = internalTimeSum/Math.max(1,internalNodes.size());
            System.out.println("Internal node time :"+internalTimeSum+", average "+internalTimeAverage);
            queue.addAll(leafNodes);
        }
        if (epsilonBound <= targetEpsilon)
            return;
        loopCleanup(newNodes, loopWatch, numNodes);
        debugHeap();
    }

    private void populateQueues(List<MARKStarNode> internalNodes, List<MARKStarNode> leafNodes, BigDecimal internalZ,
                                BigDecimal leafZ, BigDecimal[] ZSums) {
        List<MARKStarNode> leftoverLeaves = new ArrayList<>();
        int maxMinimizations = parallelism.numThreads;
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
            if (node.getLevel() < RCs.getNumPos()) {
                internalNodes.add(curNode);
                internalZ = internalZ.add(diff);
            }
            else if(shouldMinimize(node) && !correctedNode(leftoverLeaves, curNode, node)) {
                if(leafNodes.size() < maxMinimizations)
                    leafNodes.add(curNode);
                else
                    leftoverLeaves.add(curNode);
                leafZ = leafZ.add(diff);
            }

        }
        ZSums[0] = internalZ;
        ZSums[1] = leafZ;
        synchronized (this) {
            queue.addAll(leftoverLeaves);
        }
    }

    private void tightenBound() {
        System.out.println(String.format("Current overall error bound: %12.6f",epsilonBound));
        List<MARKStarNode> newNodes = new ArrayList<>();
        List<MARKStarNode> newNodesToMinimize = new ArrayList<>();
        loopPartialTime = 0;
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        double bestLower = Double.POSITIVE_INFINITY;
        if (!queue.isEmpty()) {
            bestLower = queue.peek().getConfSearchNode().getConfLowerBound();
        }
        else debugPrint("Out of conformations.");
        int numMinimizations = 0;
        int maxMinimizations = parallelism.numThreads;
        int maxNodes = 1000;
        int numNodes = 0;
        double energyThreshhold = 1000000;// -10*Math.log(((1-epsilonBound)/(1-targetEpsilon)));
        while(numMinimizations < maxMinimizations && numNodes < maxNodes &&
                !queue.isEmpty() && queue.peek().getConfSearchNode().getConfLowerBound() - bestLower < energyThreshhold ) {
            //System.out.println("Current overall error bound: "+epsilonBound);
            if(epsilonBound <= targetEpsilon)
                return;
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            double curLower = node.getConfLowerBound();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore = node.getConfLowerBound() - node.gscore;
            double confCorrection = Math.min(correctgscore, node.rigidScore) + hscore;
            if(!node.isMinimized() && node.getConfLowerBound() - confCorrection > 1e-5) {
//                debugPrint("Correcting :[" + SimpleConfSpace.formatConfRCs(node.assignments)
//                        + ":" + node.gscore + "] down to " + confCorrection);
                recordCorrection(node.getConfLowerBound(), correctgscore - node.gscore);

                node.gscore = correctgscore;
                if (confCorrection > node.rigidScore) {
                    double rigid = rigidEmat.confE(node.assignments);
                    System.out.println("Overcorrected"+SimpleConfSpace.formatConfRCs(node.assignments)+": " + confCorrection + " > " + node.rigidScore);
                    node.gscore = node.rigidScore;
                    confCorrection = node.rigidScore + hscore;
                }
                node.setBoundsFromConfLowerAndUpper(confCorrection, node.getConfUpperBound());
                curNode.markUpdated();
                queue.add(curNode);
                continue;
            }
            debugPrint("Processing Node: " + node.toString());

            //If the child is a leaf, calculate n-body minimized energies
            if (node.getLevel() < RCs.getNumPos()) {
                processPartialConfNode(newNodes, curNode, node);
            }
            else if(shouldMinimize(node) && !correctedNode(newNodes, curNode, node)) {
                synchronized (this) {
                    numMinimizations++;
                }
                processFullConfNode(newNodes, curNode, node);
                if (epsilonBound <= targetEpsilon)
                    return;
            }
            synchronized (this) {
                numNodes++;
                bestLower = Math.min(bestLower, curLower);
            }

        }

        minimizingEcalc.tasks.waitForFinish();
        internalTasks.waitForFinish();
        loopCleanup(newNodes, loopWatch, numNodes);

    }

    private void loopCleanup(List<MARKStarNode> newNodes, Stopwatch loopWatch, int numNodes) {
        synchronized (this) {
            queue.addAll(newNodes);
        }
        loopWatch.stop();
        double loopTime = loopWatch.getTimeS();
        System.out.println("Processed "+numNodes+" this loop, spawning "+newNodes.size()+" in "+loopTime+", "+stopwatch.getTime()+" so far");
        loopWatch.reset();
        loopWatch.start();
        processPreminimization(minimizingEcalc);
        System.out.println("Preminimization time : "+loopWatch.getTime(2));
        double curEpsilon = epsilonBound;
        //rootNode.updateConfBounds(new ConfIndex(RCs.getNumPos()), RCs, gscorer, hscorer);
        updateBound();
        loopWatch.stop();
        cleanupTime = loopWatch.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
        System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]",rootNode.getLowerBound(),rootNode.getUpperBound()));
    }

    private boolean correctedNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
        assert(curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            debugCorrection(node, confCorrection, oldg);
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.rigidScore);
            curNode.markUpdated();
            newNodes.add(curNode);
            return true;
        }
        return false;
    }

    private void processPartialConfNodeSerial(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
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
                    ScoreContext context = contexts.checkout();
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
                        rigiddiff=rigiddiff-node.gscore+node.rigidScore;
                        child.rigidScore = rigiddiff;

                        double confLowerBound = child.gscore + hdiff;
                        double confUpperbound = rigiddiff + maxhdiff;
                        child.computeNumConformations(RCs);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if(diff < confCorrection) {
                            /*
                            debugPrint("Correcting node " + SimpleConfSpace.formatConfRCs(child.assignments)
                                    + ":" + diff + "->" + confCorrection);
                                    */
                            recordCorrection(confLowerBound, confCorrection - diff);
                            confLowerBound = confCorrection + hdiff;
                        }
                        child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                        //progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        confRigid=confRigid-node.gscore+node.rigidScore;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if(lowerbound < confCorrection) {
                            //debugPrint("Correcting node " + SimpleConfSpace.formatConfRCs(child.assignments)
                            //        + ":" + lowerbound + "->" + confCorrection);
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection,confRigid);
                        child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                        child.gscore = child.getConfLowerBound();
                        child.rigidScore = confRigid;
                        numConfsScored++;
                        progress.reportLeafNode(child.gscore, queue.size(), epsilonBound);
                    }
                    partialTime.stop();
                    loopPartialTime+=partialTime.getTimeS();



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
        }
    }

    private void processPartialConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {
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

            internalTasks.submit(() -> {

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
                            /*
                            debugPrint("Correcting node " + SimpleConfSpace.formatConfRCs(child.assignments)
                                    + ":" + diff + "->" + confCorrection);
                                    */
                            recordCorrection(confLowerBound, confCorrection - diff);
                            confLowerBound = confCorrection + hdiff;
                        }
                        child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                        //progress.reportInternalNode(child.level, child.gscore, child.getHScore(), queue.size(), children.size(), epsilonBound);
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.rigidscorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        confRigid=confRigid-node.gscore+node.rigidScore;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if(lowerbound < confCorrection) {
                            //debugPrint("Correcting node " + SimpleConfSpace.formatConfRCs(child.assignments)
                            //        + ":" + lowerbound + "->" + confCorrection);
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection,confRigid);
                        child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                        child.gscore = child.getConfLowerBound();
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


    private void processFullConfNode(List<MARKStarNode> newNodes, MARKStarNode curNode, Node node) {

        double confCorrection = correctionMatrix.confE(node.assignments);
        if(node.getConfLowerBound() < confCorrection || node.gscore < confCorrection) {
            double oldg = node.gscore;
            node.gscore = confCorrection;
            recordCorrection(oldg, confCorrection - oldg);
            debugCorrection(node, confCorrection, oldg);
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.rigidScore);
            curNode.markUpdated();
            newNodes.add(curNode);
            return;
        }
        Stopwatch minimizationTime = new Stopwatch().start();
        leafTasks.submit(() -> {
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                node.index(context.index);

                ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound());
                ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                computeEnergyCorrection(analysis, conf, context.ecalc);

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
                double correctedEnergy = correctionMatrix.confE(node.assignments);
                recordCorrection(oldgscore, newConfLower-oldgscore);
                String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound()) + "," + (node.getConfUpperBound()) + "]";
                debugPrint(out);
                synchronized(this) {
                    numConfsEnergied++;
                }
                printMinimizationOutput(node, newConfLower, oldgscore);


            }
            return null;
        },
                // Dummy function. We're not doing anything here.
                (Node child) -> {
                    progress.reportLeafNode(node.gscore, queue.size(), epsilonBound);
                    if(!node.isMinimized())
                        newNodes.add(curNode);
                    minimizationTime.stop();
                    leafTimeAverage = Math.max(leafTimeAverage,minimizationTime.getTimeS());

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

    private void debugCorrection(Node node, double confCorrection, double oldg) {
//        debugPrint("Correcting :[" + SimpleConfSpace.formatConfRCs(node.assignments)
//                + ":" + oldg + "] down to " + confCorrection);
        if (confCorrection > node.rigidScore)
            System.err.println("Overcorrected"+SimpleConfSpace.formatConfRCs(node.assignments)+": " + confCorrection + " > " + node.rigidScore);
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
        // TODO: Replace the sortedPairwiseTerms with an ArrayList<TupE>.
        //System.out.println("Analysis:"+analysis);
        EnergyMatrix energyAnalysis = analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All);
        EnergyMatrix scoreAnalysis = analysis.breakdownScoreByPosition(minimizingEmat);
        Stopwatch correctionTime = new Stopwatch().start();
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
        double threshhold = 0.1;
        double minDifference = 0.9;//sortedPairwiseTerms.get(0).getValue()*0.7;
        double maxDiff = sortedPairwiseTerms.get(0).getValue();
        for(int i = 0; i < sortedPairwiseTerms.size(); i++)
        {
            Pair<Pair<Integer, Integer>, Double> pairEnergy = sortedPairwiseTerms.get(i);
            double pairDiff = pairEnergy.getValue();
            if(pairDiff < minDifference && maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, pairEnergy.getValue());
            int pos1 = pairEnergy.getKey().getKey();
            int pos2 = pairEnergy.getKey().getValue();
            int localMinimizations = 0;
            for(int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                computeDifference(tuple, ecalc);
                localMinimizations++;
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        ecalc.tasks.waitForFinish();
        correctionTime.stop();
        System.out.println("Correction time: "+correctionTime.getTime(2));
    }




    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if(correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        ecalc.calcEnergyAsync(tuple, (tripleEnergy)->
        {
            double lowerbound = minimizingEmat.getInternalEnergy(tuple);
            minList.set(tuple.size()-1,minList.get(tuple.size()-1)+1);
            if (tripleEnergy.energy - lowerbound > 0) {
                double correction = tripleEnergy.energy - lowerbound;
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
        int maxMinimizations = parallelism.numThreads;
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
            if (confTuple.size() > 2 && confTuple.size() < RCs.getNumPos ()){
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(minimizingEcalc, conf.toTuple());
                    return null;
                }, (econf) -> {
                });
            }
        }
        ConfIndex index = new ConfIndex(RCs.getNumPos());
        if(overlap.size() > 3 && !correctionMatrix.hasHigherOrderTermFor(overlap)
                && minimizingEmat.getInternalEnergy(overlap) != rigidEmat.getInternalEnergy(overlap)) {
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(ecalc, overlap);
                    return null;
                }, (econf) -> {
                });
        }
        minimizingEcalc.tasks.waitForFinish();
        /*
        ConfIndex index = new ConfIndex(RCs.getNumPos());
        if(overlap.size() > 3 && !correctionMatrix.hasHigherOrderTermFor(overlap)
            && minimizingEmat.getInternalEnergy(overlap) != rigidEmat.getInternalEnergy(overlap)) {
            computeTupleCorrection(ecalc, overlap);
            for (MARKStarNode conf : topConfs) {
                tasks.submit(() -> {
                    Node child = conf.getConfSearchNode();
                    child.index(index);
                    double hscore = child.getConfLowerBound() - child.gscore;
                    double confCorrection = correctionMatrix.confE(child.assignments)+hscore;
                    double confLowerBound = child.getConfLowerBound();
                    if (confLowerBound < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - child.gscore);
                        double tighterLower = confCorrection;
                        debugPrint("Correcting node " + SimpleConfSpace.formatConfRCs(child.assignments)
                                + ":" + confLowerBound + "->" + tighterLower);
                        child.setBoundsFromConfLowerAndUpper(tighterLower, child.getConfUpperBound());
                    }
                    return null;
                },(ignored)->{});
            }
        }
        */


        queue.addAll(topConfs);
    }

    private void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap) {
        if(correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        minList.set(overlap.size()-1,minList.get(overlap.size()-1)+1);
        debugPrint("Computing correction for " + overlap.stringListing() + " penalty of " + (partiallyMinimizedLower - pairwiseLower));
        progress.reportPartialMinimization(1, epsilonBound);
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

    private void debugHeap() {
        if(debug) {
            if(queue.isEmpty())
                return;
            debugPrint("A* Heap:");
            PriorityQueue<MARKStarNode> copy = new PriorityQueue<>();
            BigDecimal peek = queue.peek().getErrorBound();
            int numConfs = 0;
            while(!queue.isEmpty() && copy.size() < 10) {
                MARKStarNode node = queue.poll();
                if(node != null && (numConfs < 10 || node.getConfSearchNode().isLeaf()
                        || node.getErrorBound().multiply(new BigDecimal(1000)).compareTo(peek) >= 0)) {
                    debugPrint(node.getConfSearchNode().toString() + "," + String.format("%12.6e", node.getErrorBound()));
                    copy.add(node);
                }
                numConfs++;
            }
            queue.addAll(copy);
        }
    }

    private void updateBound() {
        double curEpsilon = epsilonBound;
        Stopwatch time = new Stopwatch().start();
        epsilonBound = rootNode.computeEpsilonErrorBounds();
        time.stop();
        System.out.println("Bound update time: "+time.getTime(2));
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
