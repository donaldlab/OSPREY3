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
import java.util.*;

public class MARKStarNode implements Comparable<MARKStarNode> {

    boolean debug = true;
    private static AStarScorer gScorer;
    private static AStarScorer rigidgScorer;
    private static AStarScorer hScorer;
    private static AStarScorer negatedHScorer;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private double errorBound = 1;
    private List<MARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode;
    public final int level;
    private static ExpFunction ef = new ExpFunction();
    private RCs RCs;


    private MARKStarNode(Node confNode) {
        confSearchNode = confNode;
        level = confSearchNode.getLevel();
        computeEpsilonErrorBounds();
        errorBound = getErrorBound();
    }

    private BigInteger getNumConfs()
    {
        return confSearchNode.numConfs;
    }


    public double computeEpsilonErrorBounds() {
        double epsilonBound = 0;
        BigDecimal lastUpper = confSearchNode.subtreeUpperBound;
        BigDecimal lastLower = confSearchNode.subtreeLowerBound;
        if(children != null && children.size() > 0) {
            BigDecimal errorUpperBound = BigDecimal.ZERO;
            BigDecimal errorLowerBound = BigDecimal.ZERO;
            for(MARKStarNode child: children) {
                double childEpsilon = child.computeEpsilonErrorBounds();
                if (debug){
                    System.out.println("Child:" + child.confSearchNode + ":" + childEpsilon);
                    System.out.println("Upper: " + setSigFigs(errorUpperBound) + "+"
                            + setSigFigs(child.confSearchNode.subtreeUpperBound) + "="
                            + setSigFigs(errorUpperBound.add(errorUpperBound)));
                    System.out.println("Lower: " + setSigFigs(errorLowerBound) + "+"
                            + setSigFigs(child.confSearchNode.subtreeLowerBound) + "="
                            + setSigFigs(errorLowerBound.add(errorLowerBound)));
                }
                errorUpperBound = errorUpperBound.add(child.confSearchNode.subtreeUpperBound);
                errorLowerBound = errorLowerBound.add(child.confSearchNode.subtreeLowerBound);
            }
            confSearchNode.subtreeUpperBound = errorUpperBound;
            confSearchNode.subtreeLowerBound = errorLowerBound;
        }
        if(confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound).compareTo(BigDecimal.ONE)<1)
        {
            return 0;
        }
        epsilonBound = confSearchNode.subtreeUpperBound.subtract(confSearchNode.subtreeLowerBound)
                .divide(confSearchNode.subtreeUpperBound,RoundingMode.HALF_UP).doubleValue();
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
        if(errorBound < epsilonBound && epsilonBound - errorBound > 0.0001) {
            System.err.println("Epsilon got bigger. Error.");
            System.err.println("UpperBound change: "+confSearchNode.subtreeUpperBound.subtract(lastUpper));
            System.err.println("LowerBound change: "+confSearchNode.subtreeLowerBound.subtract(lastLower));
        }

        errorBound = epsilonBound;
        return epsilonBound;
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

    public void printTree(String prefix, FileWriter writer)
    {
        String out = prefix+confSearchNode.confToString()+": ["+setSigFigs(confSearchNode.subtreeLowerBound)+","+setSigFigs(confSearchNode.subtreeUpperBound)+"]\n";
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
            FileWriter writer = new FileWriter(new File("ConfTreeBounds.txt"));
            printTree("",  writer);
            writer.flush();
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

    public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
        confSearchNode.updateConfLowerBound(lowerBound);
        confSearchNode.updateConfUpperBound(upperBound);
        confSearchNode.setBoundsFromConfLowerAndUpper(lowerBound, upperBound);
    }


    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }



    public static MARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix rigidEnergyMatrix,
                                        EnergyMatrix minimizingEnergyMatrix, RCs rcs,
                                        ScorerFactory gScorerFactory, ScorerFactory hscorerFactory,
                                        boolean reportProgress) {


		// make the A* scorers
		gScorer = gScorerFactory.make(minimizingEnergyMatrix);                                            // TODO: I think I want this to be minimizing
        rigidgScorer = gScorerFactory.make(rigidEnergyMatrix);
		hScorer = hscorerFactory.make(minimizingEnergyMatrix);                                            // TODO: I think I want this to be minimizing
		negatedHScorer = hscorerFactory.make(new NegatedEnergyMatrix(confSpace, rigidEnergyMatrix)); // TODO: I think I want this to be rigid

		ConfIndex confIndex = new ConfIndex(confSpace.positions.size());

		// make the root node
        Node node = new Node(confSpace.positions.size());
        Node rootNode = node;
		rootNode.index(confIndex);
		rootNode.gscore = gScorer.calc(confIndex, rcs);
        double confUpperBound = rigidgScorer.calc(confIndex,rcs)-negatedHScorer.calc(confIndex, rcs);
        double confLowerBound = rootNode.gscore+hScorer.calc(confIndex, rcs);
        rootNode.computeNumConformations(rcs);
        rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound);
		return new MARKStarNode(rootNode);
	}




    @Override
    public int compareTo(MARKStarNode other){
        return -Double.compare(this.getErrorBound(),other.getErrorBound());

    }

    public double getErrorBound() {
        if(children == null || children.size() < 1) {

            ExpFunction ef = new ExpFunction();
            Node child = confSearchNode;
            if(child.getMaxScore() > 0) {
                double diff = ef.log(ef.exp(-child.confLowerBound).subtract(ef.exp(-child.confUpperBound))).doubleValue();
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
        private BigDecimal subtreeLowerBound = null; //\hat h^ominus(f) - the lower bound on subtree contrib to partition function
        private BigDecimal subtreeUpperBound = null; //\hat h^oplus(f) - the lower bound on subtree contrib to partition function
        private double confLowerBound = -Double.MAX_VALUE;
        private double confUpperBound = Double.MAX_VALUE;
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

        public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
            updateConfLowerBound(lowerBound);
            updateConfUpperBound(upperBound);
        }


        public void updateConfLowerBound(double tighterLower) {
            if(tighterLower < confLowerBound)
                System.err.println("Updating conf lower bound of  "+confLowerBound
                        + " with "+tighterLower+", which is lower!?");
            confLowerBound = tighterLower;
            updateSubtreeUpperBound(ef.exp(-confLowerBound).multiply(new BigDecimal(getNumConformations())));
        }

        public void updateConfUpperBound(double tighterUpper) {
            if(tighterUpper > confUpperBound)
                System.err.println("Updating conf greater bound of  "+confUpperBound
                        + " with "+tighterUpper+", which is greater!?");
            confUpperBound = tighterUpper;
            updateSubtreeLowerBound(ef.exp(-confUpperBound).multiply(new BigDecimal(getNumConformations())));
        }

        public void updateSubtreeLowerBound(BigDecimal tighterLower) {
            if(subtreeLowerBound != null && subtreeLowerBound.compareTo(tighterLower)>0)
                System.err.println("Updating subtree lower bound "+setSigFigs(subtreeLowerBound)
                    + " with "+tighterLower+", which is lower!?");
            subtreeLowerBound = tighterLower;
        }

        public void updateSubtreeUpperBound(BigDecimal tighterUpper) {
            double logOriginal = 0;
            if(subtreeUpperBound != null)
                logOriginal = ef.log10(subtreeUpperBound);
            double logTighter = ef.log10(tighterUpper);
            if(subtreeUpperBound != null && subtreeUpperBound.compareTo(tighterUpper)<0)
                System.err.println("Updating subtree lower bound "+setSigFigs(subtreeUpperBound)
                        + " with "+setSigFigs(tighterUpper)+", which is lower!?");
            ExpFunction ef = new ExpFunction();
            subtreeUpperBound = tighterUpper;
        }

        public boolean isMinimized() {
            return confLowerBound == confUpperBound;
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

        public double getMaxScore() {
            return -confLowerBound;
        }

        public double getMinScore() {
            return -confUpperBound;
        }

        @Override
        public double getHScore() {
            //TODO: Scale by number of conformations
            BigDecimal d = new BigDecimal(numConfs);
            ExpFunction ef = new ExpFunction();

            return -ef.log(ef.exp(-confLowerBound).subtract(ef.exp(-confUpperBound)).multiply(d)).doubleValue()
                    -ef.log(new BigDecimal(getNumConformations())).doubleValue();
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
            out+="Energy:"+gscore+" in ["+confLowerBound+","+confUpperBound+"]->["+ setSigFigs(subtreeLowerBound) +","+ setSigFigs(subtreeUpperBound)+"]";
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
