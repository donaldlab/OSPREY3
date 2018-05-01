package edu.duke.cs.osprey.markstar.framework;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.NegatedEnergyMatrix;
import edu.duke.cs.osprey.tools.ExpFunction;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.*;

public class MARKStarNode implements Comparable<MARKStarNode> {

    private static AStarScorer gScorer;
    private static AStarScorer hScorer;
    private static AStarScorer negatedHScorer;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private BigDecimal errorUpperBound; // Note that this is actually an upper bound on pfunc of subtree
    private BigDecimal errorLowerBound;
    private double errorBound = 1;
    private List<MARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode;
    public final int level;
    private static ExpFunction ef = new ExpFunction();


    private MARKStarNode(Node confNode) {
        confSearchNode = confNode;
        errorUpperBound = ef.exp(confSearchNode.maxHScore);
        errorLowerBound = ef.exp(confSearchNode.minHScore);
        level = confSearchNode.getLevel();
        errorBound = getErrorBound();
    }

    private BigInteger getNumConfs()
    {
        return confSearchNode.numConfs;
    }


    public double computeEpsilonErrorBounds() {
        double epsilonBound = 0;
        if(children != null && children.size() > 0) {
            errorUpperBound = BigDecimal.ZERO;
            errorLowerBound = BigDecimal.ZERO;
            for(MARKStarNode child: children) {
                child.computeEpsilonErrorBounds();
                errorUpperBound = errorUpperBound.add(child.errorUpperBound);
                errorLowerBound = errorLowerBound.add(child.errorLowerBound);
            }
        }
        if(errorUpperBound.subtract(errorLowerBound).compareTo(BigDecimal.ONE)<1)
        {
            return 0;
        }
        epsilonBound = errorUpperBound.subtract(errorLowerBound).divide(errorUpperBound,RoundingMode.HALF_UP).doubleValue();
        return epsilonBound;
    }

    public BigDecimal getUpperBound(){
        return errorUpperBound;
    }

    public BigDecimal getLowerBound(){
        return errorLowerBound;
    }

    public BigDecimal setSigFigs(BigDecimal decimal, int numSigFigs)
    {
       return decimal.setScale(4-decimal.precision()+decimal.scale(),RoundingMode.HALF_UP);
    }

    public BigDecimal setSigFigs(BigDecimal decimal){
        return setSigFigs(decimal, 4);
    }

    public void printTree(String prefix, FileWriter writer)
    {
        String out = prefix+confSearchNode.confToString()+": ["+setSigFigs(errorUpperBound)+","+setSigFigs(errorLowerBound)+"]\n";
        //System.out.println(out);
        try {
            writer.write(out);
        } catch (IOException e) {
            e.printStackTrace();
        }
        if(children != null && !children.isEmpty())
           for(MARKStarNode child: children)
               child.printTree(prefix+"~+",writer);


    }

    public void printTree()
    {
        try {
            printTree("",  new FileWriter(new File("ConfTreeBounds.txt")));
        } catch (IOException e) {
            e.printStackTrace();
        }
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
                                        ScorerFactory gScorerFactory, ScorerFactory hscorerFactory,
                                        boolean reportProgress) {


		// make the A* scorers
		gScorer = gScorerFactory.make(energyMatrix);
		hScorer = hscorerFactory.make(energyMatrix);
		negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, energyMatrix));

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
		Node rootNode = new Node(confSpace.positions.size());
		rootNode.index(confIndex);
		rootNode.gscore = gScorer.calc(confIndex, rcs);
		rootNode.minHScore = hScorer.calc(confIndex, rcs);
		rootNode.maxHScore = -negatedHScorer.calc(confIndex, rcs);
        double logMax = -(rootNode.gscore+rootNode.minHScore);
        double logMin = -(rootNode.gscore+rootNode.maxHScore);
        rootNode.minHScore = logMin;
        rootNode.maxHScore = logMax;
		return new MARKStarNode(rootNode);
	}




    @Override
    public int compareTo(MARKStarNode other){
        return -Double.compare(this.errorBound,other.errorBound);

    }

    public double getErrorBound() {
        if(children == null || children.size() < 1) {

            ExpFunction ef = new ExpFunction();
            Node child = confSearchNode;
            if(child.getMaxScore() > 0) {
                double diff = ef.log(ef.exp(child.maxHScore).subtract(ef.exp(child.minHScore))).doubleValue();
                errorBound = diff;
            }
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
        public final int level;
        public BigInteger numConfs = BigInteger.ZERO;

        public Node(int size) {
            this(size,0);
        }

        public Node(int size, int level) {
            assignments = new int[size];
            Arrays.fill(assignments, Unassigned);
            this.level = level;
        }

        @Override
        public Node assign(int pos, int rc) {
            Node node = new Node(assignments.length, level+1);
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
            return minHScore;
        }

        public double getMaxScore() {
            return maxHScore;
        }

        @Override
        public double getHScore() {
            //TODO: Scale by number of conformations
            BigDecimal d = new BigDecimal(numConfs);
            ExpFunction ef = new ExpFunction();
            BigDecimal upper = ef.exp(maxHScore);
            BigDecimal lower = ef.exp(minHScore);
            BigDecimal difference = upper.subtract(lower);
            BigDecimal scaled = difference.multiply(d);
            BigDecimal log = ef.log(scaled);

            return -ef.log(ef.exp(maxHScore).subtract(ef.exp(minHScore)).multiply(d)).doubleValue();
        }

        @Override
        public void setHScore(double val) {
            throw new UnsupportedOperationException();
        }

        @Override
        public int getLevel() {
            return  level;
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

        public String confToString()
        {
            String out = "(";
            for(int i =0; i < assignments.length; i++)
            {
                out+=assignments[i]+", ";
            }
            out+=")";
            return out;
        }
        public String toString()
        {
            String out = confToString();
            out+="Conf energy:"+gscore+", max: "+maxHScore+", min: "+minHScore;
            return out;
        }

        public BigInteger getNumConformations() {
            return numConfs;
        }

        public void computeNumConformations(RCs rcs) {
            BigInteger numConfs = BigInteger.ONE;
            for (int pos=0; pos<assignments.length; pos++) {
                if (assignments[pos] == Unassigned) {
                    numConfs = numConfs.multiply(BigInteger.valueOf(rcs.getNum(pos)));
                }
            }
            this.numConfs = numConfs;
        }
    }
}
