package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;
import java.util.Map;
import java.util.concurrent.CopyOnWriteArrayList;


public class MultiSequenceSHARKStarNode implements Comparable<MultiSequenceSHARKStarNode> {


    static boolean debug = false;
    private boolean updated = true;
    /**
     * TODO: 1. Make MARKStarNodes spawn their own Node and MARKStarNode children.
     * TODO: 2. Make MARKStarNodes compute and update bounds correctly
     */

    private BigDecimal errorBound = BigDecimal.ONE;
    private double nodeEpsilon = 1;
    private MultiSequenceSHARKStarNode parent;
    private List<MultiSequenceSHARKStarNode> children; // TODO: Pick appropriate data structure
    private Node confSearchNode;
    public final int level;
    private static ExpFunction ef = new ExpFunction();
    private static BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private boolean partOfLastBound = false;

    // Information for MultiSequence SHARK* Nodes
    private Map<Sequence, MathTools.BigDecimalBounds> sequenceBounds = new HashMap<>();
    private Map<Sequence, List<MultiSequenceSHARKStarNode>> childrenMap = new HashMap<>(); // probably should override the children list


    private MultiSequenceSHARKStarNode(Node confNode, MultiSequenceSHARKStarNode parent){
        confSearchNode = confNode;
        this.level = confSearchNode.getLevel();
        this.children = new ArrayList<>();
        this.parent = parent;
        // This stuff doesn't work as well when we have multisequence stuff, moved to makeChild
        //computeEpsilonErrorBounds();
        //errorBound = getErrorBound();
    }

    /**
     * Makes a SHARKStarNode generated during flexible precomputation compatible with a new conformation space
     *
     * @param permutation   The permutation array that maps the precomputed flexible positions to the positions in the new confSpace
     *                      so, permutation[i] gives the new ConfSpace index for residue i.
     * @param size  The size of the new confSpace
     */
    public void makeNodeCompatibleWithConfSpace(int[] permutation, int size, RCs RCs){
        // first change the confSearch information

        // we copy over the new RCs based on permutation information
        Node newNode = new Node(size, this.level);
        for (int i =0; i < this.getConfSearchNode().assignments.length; i++){
            newNode.assignments[permutation[i]] = this.getConfSearchNode().assignments[i];
        }
        // Now I'm going to be hacky and just copy over the assignments
        this.getConfSearchNode().assignments = newNode.assignments;
        if (this.getConfSearchNode().pos != -1){
            this.getConfSearchNode().pos = permutation[this.getConfSearchNode().pos];
        }

        // Compute the number of conformations
        this.getConfSearchNode().computeNumConformations(RCs);
        this.updated = true;
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

    private void printBoundBreakDown(Sequence seq)
    {
        printBoundBreakDown(seq, "");
    }

    private void printBoundBreakDown(Sequence seq, String prefix)
    {
        if(level == 0) {
            System.out.println("=====================BEGIN TREE INFO==================================");
            System.out.println(prefix + confSearchNode + ": [" + setSigFigs(getSequenceBounds(seq).lower)
                    + "," + setSigFigs(getSequenceBounds(seq).upper) + "], errorBound =" + String.format("%3.3e",errorBound));
        }

        List< MultiSequenceSHARKStarNode> childrenForSequence = getChildren(seq);

        if( childrenForSequence  != null &&  childrenForSequence .size() > 0) {
            BigDecimal upper = BigDecimal.ZERO;
            BigDecimal lower = BigDecimal.ZERO;
            Collections.sort( childrenForSequence );
            prefix+="+~~";
            for(MultiSequenceSHARKStarNode child:  childrenForSequence ) {
                BigDecimal childUpper = child.getSequenceBounds(seq).upper;
                BigDecimal childLower = child.getSequenceBounds(seq).lower;
                System.out.print(prefix+child.confSearchNode+": ["+setSigFigs(childLower)
                        +","+setSigFigs(childUpper)+"], epsilon="+String.format("%3.3e",errorBound));
                System.out.print("Upper: " + setSigFigs(upper) + " + "
                        + setSigFigs(childUpper) + " = "
                        + setSigFigs(upper.add(childUpper)));
                System.out.println("Lower: " + setSigFigs(lower) + " + "
                        + setSigFigs(childLower) + " = "
                        + setSigFigs(lower.add(childLower)));
                upper = upper.add(childUpper);
                lower = lower.add(childLower);
                child.printBoundBreakDown(seq, prefix);
            }
        }
        if(level == 0) {
            System.out.println("=====================END TREE INFO==================================");
        }
    }

    public void updateSubtreeBounds(Sequence seq) {
        List<MultiSequenceSHARKStarNode> childrenForSequence = getChildren(seq);
        if(!updated)
            return;
        updated = false;
        if(childrenForSequence != null && childrenForSequence.size() > 0) {
            BigDecimal errorUpperBound = BigDecimal.ZERO;
            BigDecimal errorLowerBound = BigDecimal.ZERO;
            for(MultiSequenceSHARKStarNode child: childrenForSequence) {
                child.updateSubtreeBounds(seq);
                errorUpperBound = errorUpperBound.add(child.getSequenceBounds(seq).upper);
                errorLowerBound = errorLowerBound.add(child.getSequenceBounds(seq).lower);
            }
            getSequenceBounds(seq).upper = errorUpperBound;
            getSequenceBounds(seq).lower = errorLowerBound;
        }
    }

    private double computeEpsilon(MathTools.BigDecimalBounds bounds) {
        return bounds.upper.subtract(bounds.lower)
                .divide(bounds.upper, RoundingMode.HALF_UP).doubleValue();
    }

    public double recomputeEpsilon(Sequence seq) {
        nodeEpsilon = computeEpsilon(getSequenceBounds(seq));
        return nodeEpsilon;
    }

    public int countNodesToProcess(Sequence seq) {
        List< MultiSequenceSHARKStarNode> childrenForSeq = getChildren(seq);
        if(!updated)
            return 0;
        if(updated && (childrenForSeq == null || childrenForSeq.size() <1)) {
            return 1;
        }
        int sum = 0;
        for(MultiSequenceSHARKStarNode child: childrenForSeq) {
            sum += child.countNodesToProcess(seq);
        }
        return sum;

    }

    public double computeEpsilonErrorBounds(Sequence seq) {
        if(children == null || children.size() <1) {
            return nodeEpsilon;
        }
        if(!updated)
            return nodeEpsilon;
        double epsilonBound = 0;
        BigDecimal lastUpper = getSequenceBounds(seq).upper;
        BigDecimal lastLower = getSequenceBounds(seq).lower;
        updateSubtreeBounds(seq);
        if(getSequenceBounds(seq).upper.subtract(getSequenceBounds(seq).lower).compareTo(BigDecimal.ONE)<1)
        {
            return 0;
        }
        if(level == 0) {
            epsilonBound = computeEpsilon(getSequenceBounds(seq));
            debugChecks(lastUpper, lastLower, epsilonBound, seq);
            nodeEpsilon = epsilonBound;
            if(debug)
                printBoundBreakDown(seq);
        }
        return nodeEpsilon;
    }

    public void updateConfBounds(ConfIndex index, RCs rcs, AStarScorer gscorer, AStarScorer hScorer, Sequence seq)
    {
        List< MultiSequenceSHARKStarNode> childrenForSeq = getChildren(seq);
        if(childrenForSeq == null || childrenForSeq.size() <1) {
            confSearchNode.index(index);
            double gscore = gscorer.calc(index, rcs);
            double hscore = hScorer.calc(index, rcs);
            if(gscore+hscore > confSearchNode.getConfLowerBound())
                confSearchNode.setBoundsFromConfLowerAndUpper(gscore+hscore,confSearchNode.confUpperBound);
        }

        if(childrenForSeq != null && childrenForSeq.size() > 0) {
            for(MultiSequenceSHARKStarNode child: childrenForSeq) {
                child.updateConfBounds(index, rcs, gscorer, hScorer, seq);
            }
        }
    }

    public double updateAndReportConfBoundChange(ConfIndex index, RCs rcs, AStarScorer gscorer, AStarScorer hScorer, Sequence seq)
    {
        List<MultiSequenceSHARKStarNode> childrenForSeq = getChildren(seq);
        if(childrenForSeq  == null || childrenForSeq .size() <1) {
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
        if(childrenForSeq != null && childrenForSeq.size() > 0) {
            for(MultiSequenceSHARKStarNode child: childrenForSeq) {
                sum += child.updateAndReportConfBoundChange(index, rcs, gscorer, hScorer, seq);
            }
        }
        if(sum > 0 && level == 0)
            System.out.println("Children corrected "+sum);
        return sum;
    }
    private void debugChecks(BigDecimal lastUpper, BigDecimal lastLower, double epsilonBound, Sequence seq) {
        if (!debug)
            return;
        BigDecimal tolerance = new BigDecimal(0.00001);
        if(lastUpper != null
                && getSequenceBounds(seq).upper.subtract(lastUpper).compareTo(BigDecimal.ZERO) > 0
                && getSequenceBounds(seq).upper.subtract(lastUpper).compareTo(tolerance) > 0) {
            System.err.println("Upper bound got bigger!?");
            System.err.println("Previous: "+setSigFigs(lastUpper)+", now "+setSigFigs(getSequenceBounds(seq).upper));
            System.err.println("Increased by "+lastUpper.subtract(getSequenceBounds(seq).upper));
        }
        if(lastLower != null
                && getSequenceBounds(seq).lower.subtract(lastLower).compareTo(BigDecimal.ZERO) < 0
                && lastLower.subtract(getSequenceBounds(seq).lower).compareTo(tolerance) > 0) {
            System.err.println("Lower bound got smaller!?");
            System.err.println("Decreased by "+lastLower.subtract(getSequenceBounds(seq).lower));
        }
        if(nodeEpsilon < epsilonBound && epsilonBound - nodeEpsilon > 0.0001) {
            System.err.println("Epsilon got bigger. Error.");
            System.err.println("UpperBound change: "+getSequenceBounds(seq).upper.subtract(lastUpper));
            System.err.println("LowerBound change: "+getSequenceBounds(seq).lower.subtract(lastLower));
        }

    }

    public BigDecimal getUpperBound(Sequence seq){
        throw new UnsupportedOperationException();
    }

    public BigDecimal getLowerBound(Sequence seq){
        throw new UnsupportedOperationException();
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

    public void printTree(String prefix, FileWriter writer, SimpleConfSpace confSpace, Sequence seq)
    {
        String confString = confSearchNode.confToString();
        if(confSpace != null)
            confString = confString+"->("+confSpace.formatConfRotamersWithResidueNumbers(confSearchNode.assignments)+")";
        String out = prefix+confString+":"
                +"["+confSearchNode.confLowerBound+","+confSearchNode.confUpperBound+"]->"
                +"["+setSigFigs(getSequenceBounds(seq).lower)
                +","+setSigFigs(getSequenceBounds(seq).upper)+"]"+"\n";
        if(MathTools.isLessThan(confSearchNode.getSubtreeUpperBound(), BigDecimal.ZERO))
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
        List<MultiSequenceSHARKStarNode> children = getChildren(seq);
        if( children != null && ! children.isEmpty()) {
            Collections.sort( children, (a, b)-> -a.getSequenceBounds(seq).upper
                    .compareTo(b.getSequenceBounds(seq).upper));
            for (MultiSequenceSHARKStarNode child :  children)
                child.printTree(prefix + "~+", writer, confSpace, seq);
        }


    }

    public void printTree(String name, Sequence seq) {
        printTree(name, null, seq);
    }

    public void printTree(String name, SimpleConfSpace confspace, Sequence seq)
    {
        try {
            FileWriter writer = new FileWriter(new File(name+"ConfTreeBounds.txt"));
            printTree("",  writer, confspace, seq);
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void printTree(Sequence seq) {
        printTree("", null, null, seq);
    }

    public void printTree() {
        printTree(null);
    }

    public void index(ConfIndex confIndex) {
        confSearchNode.index(confIndex);
    }

    public Node getConfSearchNode() {
        return confSearchNode;
    }

    public MultiSequenceSHARKStarNode makeChild(Node child, Sequence seq) {
        MultiSequenceSHARKStarNode newChild = new MultiSequenceSHARKStarNode(child, this);
        newChild.computeEpsilonErrorBounds(seq);
        newChild.errorBound = getErrorBound(seq);
        getChildren(seq).add(newChild);
        children.add(newChild);
        return newChild;
    }

    public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
        confSearchNode.updateConfLowerBound(lowerBound);
        confSearchNode.updateConfUpperBound(upperBound);
        confSearchNode.setBoundsFromConfLowerAndUpper(lowerBound, upperBound);
    }

    public List<MultiSequenceSHARKStarNode> getChildren() {
        return getChildren(null);
    }

    public List<MultiSequenceSHARKStarNode> getChildren(Sequence seq) {
        if(seq == null)
            return children;
        if(!childrenMap.containsKey(seq))
            childrenMap.put(seq, new ArrayList<>());
        return childrenMap.get(seq);
    }

    public boolean isLeaf() {
        return isLeaf(null);
    }

    public boolean isLeaf(Sequence seq) {
        if(seq == null)
            return children == null || children.size() < 1;
        return getChildren(seq) == null || getChildren(seq).size() < 1;
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

    public static MultiSequenceSHARKStarNode makeRoot(SimpleConfSpace confSpace, EnergyMatrix rigidEnergyMatrix,
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
        return new MultiSequenceSHARKStarNode(rootNode, null);
    }


    @Override
    public int compareTo(MultiSequenceSHARKStarNode other){
        throw new UnsupportedOperationException("You can't compare multisequence nodes without a sequence.");
    }

    public BigDecimal getErrorBound(Sequence seq) {
        if(confSearchNode.isMinimized())
            return BigDecimal.ZERO;
        if(getChildren(seq) == null || getChildren(seq).size() < 1) {
            BigDecimal diff = getSequenceBounds(seq).upper.subtract(getSequenceBounds(seq).lower);
            return  diff.multiply(new BigDecimal(confSearchNode.minimizationRatio));
        }
        BigDecimal errorSum = BigDecimal.ZERO;
        for(MultiSequenceSHARKStarNode childNode: getChildren(seq)) {
            errorSum = errorSum.add(childNode.getErrorBound(seq));
        }
        errorBound = errorSum;
        return errorBound;
    }

    protected MathTools.BigDecimalBounds getSequenceBounds(Sequence seq) {
        if(!sequenceBounds.containsKey(seq)) {
            sequenceBounds.put(seq, new MathTools.BigDecimalBounds(confSearchNode.subtreeLowerBound,
                    confSearchNode.subtreeUpperBound));
        }
        return sequenceBounds.get(seq);
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

        public void setSubtreeUpperBound(double confLowerBound){
            this.subtreeUpperBound = computeBoundsFromEnergy(confLowerBound);
        }

        public void setSubtreeLowerBound(double confUpperBound){
            this.subtreeLowerBound = computeBoundsFromEnergy(confUpperBound);
        }

        public void setConfLowerBound(double confLowerBound){
           this.confLowerBound = confLowerBound;
           setSubtreeUpperBound(confLowerBound);
        }

        public void setConfUpperBound(double confUpperBound){
            this.confUpperBound = confUpperBound;
            setSubtreeLowerBound(confUpperBound);
        }

        public void setMinimizationRatio(double v) {
            minimizationRatio = v;
        }

        public void setBoundsFromConfLowerAndUpper(double lowerBound, double upperBound) {
            if (lowerBound - upperBound > 1e-5) {
                if (debug)
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
            if (tighterLower > confLowerBound) {
                confLowerBound = tighterLower;
                updateSubtreeUpperBound(computeBoundsFromEnergy(confLowerBound));
            }
        }

        private void updateConfUpperBound(double tighterUpper) {
            if (tighterUpper < 10 && tighterUpper - confUpperBound > 1e-5)
                System.err.println("Updating conf upper bound of  " + confUpperBound
                        + " with " + tighterUpper + ", which is greater!?");
            if (tighterUpper == Double.POSITIVE_INFINITY)
                updateSubtreeLowerBound(BigDecimal.ZERO);
            if (tighterUpper < confUpperBound) {
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
            return Math.abs(confLowerBound - confUpperBound) < 1e-5;
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
