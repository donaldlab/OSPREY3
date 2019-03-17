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

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.prototype.SimpleConf;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

public class MARKStarNode implements Comparable<MARKStarNode> {

    static boolean debug = false;
    private boolean updated = true;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private BigDecimal errorBound = BigDecimal.ONE;
    private double nodeEpsilon = 1;
    private MARKStarNode parent;
    private List<MARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode;
    public final int level;
    private static ExpFunction ef = new ExpFunction();
    private static BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private RCs RCs;
    private boolean partOfLastBound = false;


    private MARKStarNode(Node confNode, MARKStarNode markStarNode) {
        confSearchNode = confNode;
        level = confSearchNode.getLevel();
        parent = markStarNode;
        computeEpsilonErrorBounds();
        errorBound = getErrorBound();
    }

    public BigInteger getNumConfs()
    {
        return confSearchNode.numConfs;
    }

    public void markUpdated()
    {
        updated = true;
        if(level > 0)
            parent.markUpdated();
    }

    public RCTuple toTuple() {
        return new RCTuple(confSearchNode.assignments);
    }

    private void printBoundBreakDown()
    {
        printBoundBreakDown("");
    }

    private void printBoundBreakDown(String prefix)
    {
        if(level == 0) {
            System.out.println("=====================BEGIN TREE INFO==================================");
            System.out.println(prefix + confSearchNode + ": [" + setSigFigs(confSearchNode.subtreeLowerBound)
                    + "," + setSigFigs(confSearchNode.subtreeUpperBound) + "], errorBound =" + String.format("%3.3e",errorBound));
        }

        if(children != null && children.size() > 0) {
            BigDecimal upper = BigDecimal.ZERO;
            BigDecimal lower = BigDecimal.ZERO;
            Collections.sort(children);
            prefix+="+~~";
            for(MARKStarNode child: children) {
                System.out.print(prefix+child.confSearchNode+": ["+setSigFigs(child.confSearchNode.subtreeLowerBound)
                        +","+setSigFigs(child.confSearchNode.subtreeUpperBound)+"], epsilon="+String.format("%3.3e",errorBound));
                System.out.print("Upper: " + setSigFigs(upper) + " + "
                        + setSigFigs(child.confSearchNode.subtreeUpperBound) + " = "
                        + setSigFigs(upper.add(child.confSearchNode.subtreeUpperBound)));
                System.out.println("Lower: " + setSigFigs(lower) + " + "
                        + setSigFigs(child.confSearchNode.subtreeLowerBound) + " = "
                        + setSigFigs(lower.add(child.confSearchNode.subtreeLowerBound)));
                upper = upper.add(child.confSearchNode.subtreeUpperBound);
                lower = lower.add(child.confSearchNode.subtreeLowerBound);
                child.printBoundBreakDown(prefix);
            }
        }
        if(level == 0) {
            System.out.println("=====================END TREE INFO==================================");
        }
    }

    public void updateSubtreeBounds() {
        if(!updated)
            return;
        updated = false;
        if(children != null && children.size() > 0) {
            BigDecimal errorUpperBound = BigDecimal.ZERO;
            BigDecimal errorLowerBound = BigDecimal.ZERO;
            for(MARKStarNode child: children) {
                child.updateSubtreeBounds();
                errorUpperBound = errorUpperBound.add(child.confSearchNode.subtreeUpperBound);
                errorLowerBound = errorLowerBound.add(child.confSearchNode.subtreeLowerBound);
            }
            confSearchNode.subtreeUpperBound = errorUpperBound;
            confSearchNode.subtreeLowerBound = errorLowerBound;
        }
    }

    public double recomputeEpsilon() {
        nodeEpsilon = confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound)
                .divide(confSearchNode.subtreeUpperBound, RoundingMode.HALF_UP).doubleValue();
        return nodeEpsilon;
    }

    public int countNodesToProcess() {
        if(!updated)
            return 0;
        if(updated && (children == null || children.size() <1)) {
            return 1;
        }
        int sum = 0;
        for(MARKStarNode child: children) {
            sum += child.countNodesToProcess();
        }
        return sum;

    }

    public double computeEpsilonErrorBounds() {
        if(children == null || children.size() <1) {
            return nodeEpsilon;
        }
        if(!updated)
            return nodeEpsilon;
        double epsilonBound = 0;
        BigDecimal lastUpper = confSearchNode.subtreeUpperBound;
        BigDecimal lastLower = confSearchNode.subtreeLowerBound;
        updateSubtreeBounds();
        if(confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound).compareTo(BigDecimal.ONE)<1)
        {
            return 0;
        }
        if(level == 0) {
            epsilonBound = confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound)
                    .divide(confSearchNode.subtreeUpperBound, RoundingMode.HALF_UP).doubleValue();
            debugChecks(lastUpper, lastLower, epsilonBound);
            nodeEpsilon = epsilonBound;
            if(debug)
                printBoundBreakDown();
        }
        return nodeEpsilon;
    }

    public void updateConfBounds(ConfIndex index, RCs rcs, AStarScorer gscorer, AStarScorer hScorer)
    {
        if(children == null || children.size() <1) {
            confSearchNode.index(index);
            double gscore = gscorer.calc(index, rcs);
            double hscore = hScorer.calc(index, rcs);
            if(gscore+hscore > confSearchNode.getConfLowerBound())
                confSearchNode.setBoundsFromConfLowerAndUpper(gscore+hscore,confSearchNode.confUpperBound);
        }

        if(children != null && children.size() > 0) {
            for(MARKStarNode child: children) {
                child.updateConfBounds(index, rcs, gscorer, hScorer);
            }
        }
    }

    public double updateAndReportConfBoundChange(ConfIndex index, RCs rcs, AStarScorer gscorer, AStarScorer hScorer)
    {
        if(children == null || children.size() <1) {
            confSearchNode.index(index);
            double gscore = gscorer.calc(index, rcs);
            double hscore = hScorer.calc(index, rcs);
            if(gscore+hscore - confSearchNode.getConfLowerBound() > 1e-5) {
                double previousLower = confSearchNode.getConfLowerBound();
                confSearchNode.setBoundsFromConfLowerAndUpper(gscore + hscore, confSearchNode.confUpperBound);
                if(gscore+hscore < -10)
                    System.out.println("Correcting "+toTuple().stringListing()+" down to "+(gscore+hscore)+" from "+previousLower
                    +", reducing it by "+(gscore+hscore - previousLower));
                return gscore+hscore - previousLower;
            }
        }
        double sum = 0;
        if(children != null && children.size() > 0) {
            for(MARKStarNode child: children) {
                sum += child.updateAndReportConfBoundChange(index, rcs, gscorer, hScorer);
            }
        }
        if(sum > 0 && level == 0)
            System.out.println("Children corrected "+sum);
        return sum;
    }
    private void debugChecks(BigDecimal lastUpper, BigDecimal lastLower, double epsilonBound) {
        if (!debug)
            return;
        BigDecimal tolerance = new BigDecimal(0.00001);
        if(lastUpper != null
                && confSearchNode.subtreeUpperBound.subtract(lastUpper).compareTo(BigDecimal.ZERO) > 0
                && confSearchNode.subtreeUpperBound.subtract(lastUpper).compareTo(tolerance) > 0) {
            System.err.println("Upper bound got bigger!?");
            System.err.println("Previous: "+setSigFigs(lastUpper)+", now "+setSigFigs(confSearchNode.subtreeUpperBound));
            System.err.println("Increased by "+lastUpper.subtract(confSearchNode.subtreeUpperBound));
        }
        if(lastLower != null
                && confSearchNode.subtreeLowerBound.subtract(lastLower).compareTo(BigDecimal.ZERO) < 0
                && lastLower.subtract(confSearchNode.subtreeLowerBound).compareTo(tolerance) > 0) {
            System.err.println("Lower bound got smaller!?");
            System.err.println("Decreased by "+lastLower.subtract(confSearchNode.subtreeLowerBound));
        }
        if(nodeEpsilon < epsilonBound && epsilonBound - nodeEpsilon > 0.0001) {
            System.err.println("Epsilon got bigger. Error.");
            System.err.println("UpperBound change: "+confSearchNode.subtreeUpperBound.subtract(lastUpper));
            System.err.println("LowerBound change: "+confSearchNode.subtreeLowerBound.subtract(lastLower));
        }

    }


    public BigDecimal getUpperBound(){
        return confSearchNode.subtreeUpperBound;
    }

    public BigDecimal getLowerBound(){
        return confSearchNode.subtreeLowerBound;
    }

    public static BigDecimal setSigFigs(BigDecimal decimal, int numSigFigs) {
       return decimal.setScale(4-decimal.precision()+decimal.scale(),RoundingMode.HALF_UP);
    }

    public static BigDecimal setSigFigs(BigDecimal decimal){
        return setSigFigs(decimal, 4);
    }


    public void printTree(String prefix, FileWriter writer, SimpleConfSpace confSpace)
    {
        String confString = confSearchNode.confToString();
        if(confSpace != null)
            confString = confString+"->("+confSpace.formatConfRotamersWithResidueNumbers(confSearchNode.assignments)+")";
        String out = prefix+confString+":"
                +"["+confSearchNode.confLowerBound+","+confSearchNode.confUpperBound+"]->"
                +"["+setSigFigs(confSearchNode.subtreeLowerBound)
                +","+setSigFigs(confSearchNode.subtreeUpperBound)+"]"+"\n";
        if(MathTools.isLessThan(confSearchNode.getSubtreeUpperBound(), BigDecimal.ONE))
            return;
        if(writer != null) {
            try {
                writer.write(out);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else
            System.out.print(out);
        if(children != null && !children.isEmpty()) {
            Collections.sort(children, (a,b)-> -a.confSearchNode.subtreeUpperBound
                    .compareTo(b.confSearchNode.subtreeUpperBound));
            for (MARKStarNode child : children)
                child.printTree(prefix + "~+", writer, confSpace);
        }


    }

    public void printTree(String name) {
        printTree(name, null);
    }

    public void printTree(String name, SimpleConfSpace confspace)
    {
        try {
            FileWriter writer = new FileWriter(new File(name+"ConfTreeBounds.txt"));
            printTree("",  writer, confspace);
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void printTree() {
        printTree("", null, null);
    }

    public void index(ConfIndex confIndex) {
        confSearchNode.index(confIndex);
    }

    public Node getConfSearchNode() {
        return confSearchNode;
    }

    public MARKStarNode makeChild(Node child) {
        MARKStarNode newChild = new MARKStarNode(child, this);
        if(children == null)
            children = new CopyOnWriteArrayList<>();
        children.add(newChild);
        return newChild;
    }

    public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
        confSearchNode.updateConfLowerBound(lowerBound);
        confSearchNode.updateConfUpperBound(upperBound);
        confSearchNode.setBoundsFromConfLowerAndUpper(lowerBound, upperBound);
    }

    public List<? extends MARKStarNode> getChildren() {
        return children;
    }

    public static enum Type {
        internal,
        boundedLeaf,
        minimizedLeaf
    }
    public Type type(RCs rcs) {
        if(level < rcs.getNumPos())
            return  Type.internal;
        if(!confSearchNode.isMinimized())
            return Type.boundedLeaf;
        return Type.minimizedLeaf;
    }





    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }



    public static MARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix rigidEnergyMatrix,
                                        EnergyMatrix minimizingEnergyMatrix, RCs rcs,
                                        AStarScorer gScorer, AStarScorer hScorer,
                                        AStarScorer rigidgScorer, AStarScorer negatedHScorer,
                                        boolean reportProgress) {


		// make the A* scorers

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
        Node node = new Node(confSpace.positions.size());
        Node rootNode = node;
		rootNode.index(confIndex);
		rootNode.gscore = gScorer.calc(confIndex, rcs);
		rootNode.rigidScore = rigidgScorer.calc(confIndex,rcs);
        double confUpperBound = rigidgScorer.calc(confIndex,rcs)-negatedHScorer.calc(confIndex, rcs);
        double confLowerBound = rootNode.gscore+hScorer.calc(confIndex, rcs);
        rootNode.computeNumConformations(rcs);
        rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound);
		return new MARKStarNode(rootNode, null);
	}


    @Override
    public int compareTo(MARKStarNode other){
        return -getErrorBound().compareTo(other.getErrorBound());
    }

    public BigDecimal getErrorBound() {
        if(confSearchNode.isMinimized())
            return BigDecimal.ZERO;
        if(children == null || children.size() < 1) {
            BigDecimal diff = confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound);
            return  diff.multiply(new BigDecimal(confSearchNode.minimizationRatio));
        }
        BigDecimal errorSum = BigDecimal.ZERO;
        for(MARKStarNode childNode: children) {
            errorSum = errorSum.add(childNode.getErrorBound());
        }
        errorBound = errorSum;
        return errorBound;
    }


    public static class Node implements ConfAStarNode {

        private static int Unassigned = -1;
        public double gscore = Double.NaN;
        public double rigidScore = Double.NaN;
        private BigDecimal subtreeLowerBound = BigDecimal.ZERO; //\hat h^ominus(f) - the lower bound on subtree contrib to partition function
        private BigDecimal subtreeUpperBound = null; //\hat h^oplus(f) - the lower bound on subtree contrib to partition function
        private double confLowerBound = -Double.MAX_VALUE;
        private double confUpperBound = Double.MAX_VALUE;
        public int[] assignments;
        public int pos = Unassigned;
        public int rc = Unassigned;
        public final int level;
        public BigInteger numConfs = BigInteger.ZERO;
        private double minimizationRatio = 1;

        public Node(int size) {
            this(size, 0);
        }

        public Node(int size, int level) {
            assignments = new int[size];
            Arrays.fill(assignments, Unassigned);
            this.level = level;
        }

        public void setMinimizationRatio(double v)
        {
            minimizationRatio = v;
        }

        public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
            if (lowerBound - upperBound > 1e-5){
                if(debug)
                    System.out.println("Incorrect conf bounds set.");
                double temp = lowerBound;
                lowerBound = upperBound;
                upperBound = temp;
                lowerBound = Math.min(0, lowerBound);
                upperBound = Math.max(lowerBound, upperBound);
            }
            updateConfLowerBound(lowerBound);
            updateConfUpperBound(upperBound);
        }


        private void updateConfLowerBound(double tighterLower) {
            if (tighterLower < 10 && tighterLower - confLowerBound < -1e-5)
                System.err.println("Updating conf lower bound of  " + confLowerBound
                        + " with " + tighterLower + ", which is lower!?");
            if(tighterLower > confLowerBound) {
                confLowerBound = tighterLower;
                updateSubtreeUpperBound(computeBoundsFromEnergy(confLowerBound));
            }
        }

        private void updateConfUpperBound(double tighterUpper) {
            if (tighterUpper < 10 && tighterUpper - confUpperBound > 1e-5)
                System.err.println("Updating conf upper bound of  " + confUpperBound
                        + " with " + tighterUpper + ", which is greater!?");
            if(tighterUpper == Double.POSITIVE_INFINITY)
                updateSubtreeLowerBound(BigDecimal.ZERO);
            if(tighterUpper < confUpperBound) {
                confUpperBound = tighterUpper;
                updateSubtreeLowerBound(computeBoundsFromEnergy(confUpperBound));
            }
        }

        private BigDecimal computeBoundsFromEnergy(double energy) {
            return bc.calc(energy).multiply(new BigDecimal(getNumConformations()));
        }

        private void updateSubtreeLowerBound(BigDecimal tighterLower) {
            if (subtreeLowerBound != null && subtreeLowerBound.compareTo(tighterLower) > 0)
                System.err.println("Updating subtree lower bound " + setSigFigs(subtreeLowerBound)
                        + " with " + tighterLower + ", which is lower!?");
            subtreeLowerBound = tighterLower;
        }

        private void updateSubtreeUpperBound(BigDecimal tighterUpper) {
            if (subtreeUpperBound != null && subtreeUpperBound.compareTo(tighterUpper) < 0)
                System.err.println("Updating subtree upper bound " + setSigFigs(subtreeUpperBound)
                        + " with " + setSigFigs(tighterUpper) + ", which is greater!?");
            subtreeUpperBound = tighterUpper;
        }

        public boolean isMinimized() {
            return confLowerBound == confUpperBound;
        }

        public boolean isLeaf() {
            return level == assignments.length;
        }

        @Override
        public Node assign(int pos, int rc) {
            Node node = new Node(assignments.length, level + 1);
            node.pos = pos;
            node.rc = rc;
            System.arraycopy(assignments, 0, node.assignments, 0, assignments.length);
            node.assignments[pos] = rc;
            return node;
        }

        @Override
        public double getGScore() {
            return gscore;
        }

        @Override
        public void setGScore(double val) {
            gscore = val;
        }

        public double getMaxScore() {
            return -confLowerBound;
        }

        public double getMinScore() {
            return -confUpperBound;
        }

        @Override
        public double getHScore() {
            return -subtreeUpperBound.subtract(subtreeLowerBound).doubleValue();
        }

        @Override
        public void setHScore(double val) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int getLevel() {
            return level;
        }

        @Override
        public void getConf(int[] conf) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void index(ConfIndex confIndex) {
            confIndex.numDefined = 0;
            confIndex.numUndefined = 0;
            for (int pos = 0; pos < assignments.length; pos++) {
                int rc = assignments[pos];
                if (rc == -1) {
                    confIndex.undefinedPos[confIndex.numUndefined] = pos;
                    confIndex.numUndefined++;
                } else {
                    confIndex.definedPos[confIndex.numDefined] = pos;
                    confIndex.definedRCs[confIndex.numDefined] = assignments[pos];
                    confIndex.numDefined++;
                }
            }
            confIndex.node = this;
        }

        public String confToString() {
            String out = "(";
            for (int i = 0; i < assignments.length; i++) {
                out += assignments[i] + ", ";
            }
            out += ")";
            return out;
        }

        public String toString() {
            String out = confToString();
            out += "Energy:" + String.format("%4.2f", gscore) + "*" + numConfs;
            if (!isMinimized())
                out += " in [" + String.format("%4.4e,%4.4e", confLowerBound, confUpperBound) + "]->[" + setSigFigs(subtreeLowerBound) + "," + setSigFigs(subtreeUpperBound) + "]";
            else
                out += " (minimized) -> " + setSigFigs(subtreeLowerBound);
            return out;
        }

        public BigInteger getNumConformations() {
            return numConfs;
        }

        public void computeNumConformations(RCs rcs) {
            BigInteger numConfs = BigInteger.ONE;
            this.numConfs = numConfs;
            if(rcs.getNumPos() == assignments.length) {
                boolean fullyAssigned = true;
                for (int pos = 0; pos < assignments.length; pos++) {
                    if(assignments[pos] == Unassigned)
                        fullyAssigned = false;
                }
                if(fullyAssigned)
                    return;
            }

            for (int pos = 0; pos < assignments.length; pos++) {
                if (assignments[pos] == Unassigned) {
                    numConfs = numConfs.multiply(BigInteger.valueOf(rcs.getNum(pos)));
                }
            }
            this.numConfs = numConfs;
            assert(this.numConfs.compareTo(BigInteger.ZERO) > 0);
        }

        public double getConfLowerBound() {
            return confLowerBound;
        }

        public double getConfUpperBound() {
            return confUpperBound;
        }

        public BigDecimal getSubtreeLowerBound() {
            return subtreeLowerBound;
        }

        public BigDecimal getSubtreeUpperBound() {
            return subtreeUpperBound;
        }
    }

}
