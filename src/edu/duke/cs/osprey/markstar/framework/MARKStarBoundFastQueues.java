package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.MSKStar;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.math.BigDecimal;
import java.util.*;

public class MARKStarBoundFastQueues extends MARKStarBound {

    public String stateName = String.format("%4f",Math.random());
    private Queue<MARKStarNode> leafQueue;
    private Queue<MARKStarNode> internalQueue;



    public static MARKStarBoundFastQueues makeFromConfSpaceInfo(SimpleConfSpace confSpace, ConfEnergyCalculator minimizingConfEcalc, EnergyMatrix ematRigid, EnergyMatrix ematMinimized, RCs rcs) {
        return new MARKStarBoundFastQueues(confSpace, ematRigid, ematMinimized, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
    }

    public static MARKStarBoundFastQueues makeFromConfSpaceInfo(BBKStar.ConfSpaceInfo info, ConfEnergyCalculator minimizingConfEcalc, RCs rcs) {
        return new MARKStarBoundFastQueues(info.confSpace, info.rigidEmat, info.minimizingEmat, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
    }

    public static MARKStarBoundFastQueues makeFromMSKState(MSKStar.State info, ConfEnergyCalculator minimizingConfEcalc, RCs rcs) {
        return new MARKStarBoundFastQueues(info.confSpace, null, null, minimizingConfEcalc, rcs, minimizingConfEcalc.ecalc.parallelism);
    }

    public MARKStarBoundFastQueues(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                   ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        super(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        this.leafQueue = new PriorityQueue<>();
        this.internalQueue = new PriorityQueue<>();
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
            if(maxMinimizations < parallelism.numThreads)
                maxMinimizations++;
            internalQueue.addAll(internalNodes);
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
            numInternalNodesProcessed+=internalNodes.size();
            leafQueue.addAll(leafNodes);
        }
        if (epsilonBound <= targetEpsilon)
            return;
        loopCleanup(newNodes, loopWatch, numNodes);
    }

    @Override
    protected void populateQueues(Queue<MARKStarNode> queue, List<MARKStarNode> internalNodes, List<MARKStarNode> leafNodes, BigDecimal internalZ,
                                BigDecimal leafZ, BigDecimal[] ZSums) {
        List<MARKStarNode> leftoverLeaves = new ArrayList<>();
        int maxNodes = 1000;
        if(leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int)Math.floor(0.1*leafTimeAverage/internalTimeAverage));
        while(!queue.isEmpty() && (internalQueue.size() < maxNodes || leafQueue.size() < maxMinimizations)){
            MARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(super.RCs.getNumPos());
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

            if (node.getLevel() < RCs.getNumPos()) {
                internalQueue.add(curNode);
            }
            else if(shouldMinimize(node) && !correctedNode(leftoverLeaves, curNode, node)) {
                leafQueue.add(curNode);
            }

        }

        ZSums[0] = fillListFromQueue(internalNodes, internalQueue, maxNodes);
        ZSums[1] = fillListFromQueue(leafNodes, leafQueue, maxMinimizations);
        queue.addAll(leftoverLeaves);
    }

    private BigDecimal fillListFromQueue(List<MARKStarNode> list, Queue<MARKStarNode> queue, int max) {
        BigDecimal sum = BigDecimal.ZERO;
        List<MARKStarNode> leftovers = new ArrayList<>();
        while(!queue.isEmpty() && list.size() < max) {
            MARKStarNode curNode = queue.poll();
            if(correctedNode(leftovers, curNode, curNode.getConfSearchNode())) {
                continue;
            }
            BigDecimal diff = curNode.getUpperBound().subtract(curNode.getLowerBound());
            sum = sum.add(diff);
            list.add(curNode);
        }
        queue.addAll(leftovers);
        return sum;
    }


}
