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
                System.out.println("Correction from "+correctionMatrix.sourceECalc+":"+node.gscore+"->"+correctgscore);
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
