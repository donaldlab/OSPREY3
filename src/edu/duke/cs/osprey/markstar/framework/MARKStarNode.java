package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;

public class MARKStarNode implements Comparable<MARKStarNode> {

    private static AStarScorer gScorer;
    private static AStarScorer hScorer;
    private static AStarScorer negatedHScorer;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private double errorUpperBound;
    private double errorLowerBound;
    private double errorBound;
    private List<MARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode;


    private MARKStarNode(Node confNode) {
        confSearchNode = confNode;
        errorUpperBound = confSearchNode.maxHScore;
        errorLowerBound = confSearchNode.minHScore;
        computeErrorBounds();
    }

    public void scoreNode(Node confNode, ConfIndex confIndex, RCs rcs) {
        confNode.gscore = gScorer.calc(confIndex, rcs);
        confNode.minHScore = hScorer.calc(confIndex, rcs);
        confNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);

    }


    private double computeErrorBounds() {
        errorBound = 0;
        System.out.println("Upper and Lower bounds : "+errorUpperBound+", "+errorLowerBound);
        if(children == null || children.size() < 1)
            errorBound = errorUpperBound-errorLowerBound;
        else {
            for(MARKStarNode child: children)
            {
                errorBound+= child.computeErrorBounds();
            }
            errorBound /= children.size();
        }
        return errorBound;
    }

    public void index(ConfIndex confIndex) {
        confSearchNode.index(confIndex);
    }

    public Node getConfSearchNode() {
        return confSearchNode;
    }

    public MARKStarNode makeChild(Node child) {
        MARKStarNode newChild = new MARKStarNode(child);
        if(children == null)
            children = new ArrayList<>();
        children.add(newChild);
        return newChild;
    }


    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }



    public static MARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix energyMatrix, RCs rcs,
                                        ScorerFactory gscorerFactory, ScorerFactory hscorerFactory,
                                        boolean reportProgress) {


		// make the A* scorers
		gScorer = gscorerFactory.make(energyMatrix);
		hScorer = hscorerFactory.make(energyMatrix);
		negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, energyMatrix));

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
		Node rootNode = new Node(confSpace.positions.size());
		rootNode.index(confIndex);
		rootNode.gscore = gScorer.calc(confIndex, rcs);
		rootNode.minHScore = hScorer.calc(confIndex, rcs);
		rootNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);
		return new MARKStarNode(rootNode);
	}




    @Override
    public int compareTo(MARKStarNode other){
        return Double.compare(this.errorBound,other.errorBound);

    }


    public double getErrorBound() {
        if(children == null || children.size() < 1) {
            errorBound = -confSearchNode.getHScore();
            return errorBound;
        }
        double errorSum = 0;
        for(MARKStarNode childNode: children) {
            errorSum += childNode.getErrorBound();
        }
        errorBound = errorSum;
        return errorBound;
    }


    public static class Node implements ConfAStarNode {

        private static int Unassigned = -1;
        public double gscore = Double.NaN;
        public double minHScore = Double.NaN;
        public double maxHScore = Double.NaN;
        public int[] assignments;
        public int pos = Unassigned;
        public int rc = Unassigned;

        public Node(int size) {
            assignments = new int[size];
            Arrays.fill(assignments, Unassigned);
        }

        @Override
        public Node assign(int pos, int rc) {
            Node node = new Node(assignments.length);
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

        public double getMinScore() {
            return gscore + minHScore;
        }

        public double getMaxScore() {
            return gscore + maxHScore;
        }

        @Override
        public double getHScore() {
            // We want it to be as small as possible since our A* implementation finds the min
            ExpFunction ef = new ExpFunction();
            BigDecimal upperBound = ef.exp(-getMinScore());
            System.out.println("g:"+gscore+", max: "+maxHScore+", min: "+minHScore);
            System.out.println("Upper bound: "+upperBound.setScale(4,RoundingMode.FLOOR).toEngineeringString());
            double ErrorBound = 0;
            if(upperBound.doubleValue() > 0)
                ErrorBound = upperBound.subtract(ef.exp(-getMaxScore())).divide(upperBound).doubleValue();

            return -ErrorBound;
        }

        @Override
        public void setHScore(double val) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int getLevel() {
            throw new UnsupportedOperationException();
        }

        @Override
        public void getConf(int[] conf) {
            throw new UnsupportedOperationException();
        }

        @Override
        public void index(ConfIndex confIndex) {
            confIndex.numDefined = 0;
            confIndex.numUndefined = 0;
            for (int pos=0; pos<assignments.length; pos++) {
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

        public BigInteger getNumConformations(RCs rcs) {
            BigInteger numConfs = BigInteger.ONE;
            for (int pos=0; pos<assignments.length; pos++) {
                if (assignments[pos] == Unassigned) {
                    numConfs = numConfs.multiply(BigInteger.valueOf(rcs.getNum(pos)));
                }
            }
            return numConfs;
        }
    }
}
