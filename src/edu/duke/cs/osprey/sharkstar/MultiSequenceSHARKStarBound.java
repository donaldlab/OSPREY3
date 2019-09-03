package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseRigidGScorer;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.StaticBiggestLowerboundDifferenceOrder;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.printTree;

public class MultiSequenceSHARKStarBound implements PartitionFunction {

    private Sequence precomputedSequence;
    protected double targetEpsilon = 1;
    public boolean debug = true;
    public boolean profileOutput = false;
    private PartitionFunction.Status status = null;

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
    public String stateName = String.format("%4f", Math.random());
    private int numPartialMinimizations;
    private ArrayList<Integer> minList;
    protected double internalTimeAverage;
    protected double leafTimeAverage;
    private double cleanupTime;
    private boolean nonZeroLower;
    protected static TaskExecutor loopTasks;


    // We keep track of the root node for computing our K* bounds
    public MultiSequenceSHARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    private ConfIndex confIndex;
    public StaticBiggestLowerboundDifferenceOrder order;
    public final AStarPruner pruner;
    // TODO: Implement new AStarPruner for MARK*?
    protected RCs fullRCs;
    protected Parallelism parallelism;
    private ObjectPool<ScoreContext> contexts;
    private ScorerFactory gscorerFactory;
    private ScorerFactory rigidgscorerFactory;
    private ScorerFactory hscorerFactory;
    private ScorerFactory nhscorerFactory;

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

    private MultiSequenceSHARKStarBound precomputedPfunc;
    public MultiSequenceSHARKStarNode precomputedRootNode;
    private SimpleConfSpace confSpace;

    private BigDecimal precomputedUpperBound;
    private BigDecimal precomputedLowerBound;

    private List<MultiSequenceSHARKStarNode> precomputedFringe = new ArrayList<>();

    private static final int[] debugConf = new int[]{};

    /**
     * Constructor to make a default SHARKStarBound Class
     *
     * @param confSpace           the partition function conformation space
     * @param rigidEmat           the rigid pairwise energy matrix
     * @param minimizingEmat      the parwise-minimized energy matrix
     * @param minimizingConfEcalc the energy calculator to calculate minimized conf energies
     * @param rcs                 information on possible rotamers at all design positions
     * @param parallelism         information for threading
     */
    public MultiSequenceSHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                       ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        this.minimizingEcalc = minimizingConfEcalc;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);
        rigidgscorerFactory = (emats) -> new PairwiseRigidGScorer(emats);


        hscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, rcs);
        nhscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, rcs, true);
        //hscorerFactory = (emats) -> new TraditionalPairwiseHScorer(emats, rcs);
        //nhscorerFactory = (emats) -> new TraditionalPairwiseHScorer(new NegatedEnergyMatrix(confSpace, rigidEmat), rcs);

        // No precomputed sequence means the "precomputed" sequence is empty
        this.precomputedSequence = confSpace.makeUnassignedSequence();
        confIndex = new ConfIndex(rcs.getNumPos());

        Node rootConfNode = new Node(confSpace.positions.size(), 0, new MathTools.DoubleBounds());
        rootConfNode.index(confIndex);
        double partialConfLowerbound = gscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        double partialConfUpperBound = rigidgscorerFactory.make(rigidEmat).calc(confIndex, rcs);
        rootConfNode.computeNumConformations(rcs);
        rootConfNode.setBoundsFromConfLowerAndUpper(partialConfLowerbound, partialConfUpperBound);

        this.rootNode = MultiSequenceSHARKStarNode.makeRoot(rootConfNode, confSpace);
        double confLowerBound = partialConfLowerbound + hscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        double confUpperBound = partialConfLowerbound - nhscorerFactory.make(rigidEmat).calc(confIndex, rcs);
        rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound, precomputedSequence);
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.fullRCs = rcs;
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat));
        this.pruner = null;

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.partialConfLowerBoundScorer = gscorerFactory.make(minimizingEmat);
            context.lowerBoundScorer = hscorerFactory.make(minimizingEmat);
            context.partialConfUpperBoundScorer = rigidgscorerFactory.make(rigidEmat);
            /** These scoreres should match the scorers in the SHARKStarNode root - they perform the same calculations**/
            context.upperBoundScorer = nhscorerFactory.make(rigidEmat); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;
            return context;
        });

        progress = new MARKStarProgress(fullRCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        setParallelism(parallelism);

        // Recording pfunc starting bounds
        this.startLowerBound = rootNode.getLowerBound(precomputedSequence);
        this.startUpperBound = rootNode.getUpperBound(precomputedSequence);
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(), 0));
        this.confSpace = confSpace;
    }

    public MultiSequenceSHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                       ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism,
                                       MultiSequenceSHARKStarBound precomputedFlex) {

        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);

		/*
		Now we need to do a couple things.
		For now, let's assume we are working with a single sequence

		 */

        processPrecomputedFlex(precomputedFlex);


		/*
		TODO: Go through the tree to make sure that the assignments are compatible with the new confspace

		TODO: Make sure all of the energy matrices (including the correction emats) are compatible with the new confspace

		TODO: Do something with score / object contexts?

		TODO: Go through and update bounds

		TODO: Populate queue
		 */
    }

    private void processPrecomputedFlex(MultiSequenceSHARKStarBound precomputedFlex) {
        precomputedPfunc = precomputedFlex;
        precomputedRootNode = precomputedFlex.rootNode;
        this.precomputedSequence = precomputedFlex.confSpace.makeWildTypeSequence();
        precomputedUpperBound = precomputedRootNode.getUpperBound(precomputedSequence);
        precomputedLowerBound = precomputedRootNode.getLowerBound(precomputedSequence);
        updatePrecomputedConfTree();

        // Fix order issues
        ConfIndex rootIndex = new ConfIndex(fullRCs.getNumPos());
        this.rootNode.getConfSearchNode().index(rootIndex);
        this.order.updateForPrecomputedOrder(precomputedFlex.order, rootIndex, this.fullRCs, genConfSpaceMapping());
    }

    /**
     * Returns a wrapped pointer to this class, so that BBK* and MSK* can pretend they have single-sequence
     * partition functions.
     */
    public PartitionFunction getPartitionFunctionForSequence(Sequence seq) {
        SingleSequenceSHARKStarBound newBound = new SingleSequenceSHARKStarBound(seq, this);
        newBound.init(null, null, targetEpsilon);
        System.out.println("Full RCs: "+fullRCs);
        System.out.println("Sequence RCs: "+newBound.seqRCs);
        computeFringeForSequence(newBound, this.rootNode);
        newBound.updateBound();
        rootNode.updateSubtreeBounds(seq);
        //printTree(seq, this.rootNode);
        return newBound;
    }

    /**
     * Returns the partition function lower bound for a particular sequence
     * <p>
     * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
     *
     * @param seq Sequence for which to get pfunc lower bound
     * @return BigDecimal pfunc lower bound
     */
    public BigDecimal getLowerBound(Sequence seq) {
        throw new UnsupportedOperationException("getLowerBound(seq) is not yet implemented");
    }

    /**
     * Returns the partition function upper bound for a particular sequence
     * <p>
     * Note that SHARKStarBound will eventually contain a multi-sequence confTree, although this isn't currently the case
     *
     * @param seq Sequence for which to get pfunc upper bound
     * @return BigDecimal pfunc upper bound
     */
    public BigDecimal getUpperBound(Sequence seq) {
        throw new UnsupportedOperationException("getUpperBound(seq) is not yet implemented");
    }


    /**
     * Returns the partition function lower bound for the precomputed confspace
     *
     * @return BigDecimal precomputed pfunc lower bound
     */
    public BigDecimal getPrecomputedLowerBound() {
        return precomputedLowerBound;
    }

    /**
     * Returns the partition function upper bound for the precomputed confTree
     *
     * @return BigDecimal precomputed pfunc upper bound
     */
    public BigDecimal getPrecomputedUpperBound() {
        return precomputedUpperBound;
    }

    /**
     * Makes the current confTree consistent with the current confSpace
     * <p>
     * When we precompute flexible residues, we will have a tree that is for a flexible confspace.
     * However, when we want to compute for mutable residues, we need to extend the length of assignments in our tree
     */
    public void updatePrecomputedConfTree() {
        int[] permutationArray = genConfSpaceMapping();
        updatePrecomputedNode(precomputedRootNode, permutationArray, this.confSpace.getNumPos());
        this.rootNode = precomputedRootNode;
        //System.out.println("The precomputed root node is " + precomputedRootNode.toTuple());
        //System.out.println("\n###############\nFull root upper: "+rootNode.getUpperBound()+" lower: "+rootNode.getLowerBound());

		/*
		Here's the plan: use the permutation matrix to map assignements onto the new tree.
		The g scores should map fine
		the h scores may need to be updated, which is kind of awkward I suppose.
		But, any full conformations should be minimized and should be added to the MAE energy matrix
		likewise for any partial minimizations
		 */

        //TODO: update partial minimizations and full conformations?
    }

    private void updatePrecomputedNode(MultiSequenceSHARKStarNode node, int[] permutation, int size) {
        if (node.getChildren(precomputedSequence) != null) {
            for (MultiSequenceSHARKStarNode child : node.getChildren(precomputedSequence)) {
                updatePrecomputedNode(child, permutation, size);
            }
        }
        node.makeNodeCompatibleWithConfSpace(permutation, size, this.fullRCs);
        node.setNewConfSpace(confSpace);
    }

    private void initForSeq(List<MultiSequenceSHARKStarNode> nodes, SingleSequenceSHARKStarBound bound) {
        for(MultiSequenceSHARKStarNode node:nodes) {
            Node confNode = node.getConfSearchNode();
            RCs rcs = bound.seqRCs;
            confNode.index(confIndex);

            double confLowerBound = confNode.getPartialConfLowerBound() + hscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
            double confUpperBound = confNode.getPartialConfLowerBound() - nhscorerFactory.make(rigidEmat).calc(confIndex, rcs);
            node.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound, bound.sequence);
            initForSeq(node.getChildren(bound.sequence), bound);
        }

    }

    /**
     * Add the newly updated nodes to the queue so that we don't redo any work
     */
    private void addPrecomputedFringeToQueues(SingleSequenceSHARKStarBound bound) {
        processPrecomputedFringe(this.rootNode, bound);
        initForSeq(precomputedFringe, bound);
        bound.fringeNodes.addAll(precomputedFringe);

    }

    private void computeFringeForSequence(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarNode curNode) {
        Node confNode = curNode.getConfSearchNode();
        RCs rcs = bound.seqRCs;
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            confNode.index(context.index);

            RCs test = new RCs(minimizingEcalc.confSpace);
            double confLowerBound = confNode.getPartialConfLowerBound() + context.lowerBoundScorer.calc(context.index, rcs);
            double confUpperBound = confNode.getPartialConfUpperBound() + context.upperBoundScorer.calc(context.index, rcs);
            curNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound, bound.sequence);
        }
        if(curNode.getChildren(bound.sequence).size() < 1)
            bound.fringeNodes.add(curNode);
        else
            for(MultiSequenceSHARKStarNode child: curNode.getChildren(bound.sequence))
                computeFringeForSequence(bound, child);
    }

    private void processPrecomputedFringe(MultiSequenceSHARKStarNode root, SingleSequenceSHARKStarBound bound) {
        if (root.isLeaf()) {
            // Compute correct hscores
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                // index the node
                root.index(context.index);
                double hscore = context.lowerBoundScorer.calc(context.index, bound.seqRCs);
                double maxhscore = context.upperBoundScorer.calc(context.index, bound.seqRCs);
                Node confNode = root.getConfSearchNode();
                double confLowerBound = confNode.getPartialConfLowerBound() + hscore;
                double confUpperBound = confNode.getPartialConfUpperBound() + maxhscore;
                root.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound, bound.sequence);
            }
            // add node to queue
            precomputedFringe.add(root);
        } else {
            for (MultiSequenceSHARKStarNode node : root.getChildren(bound.sequence)) {
                processPrecomputedFringe(node, bound);
            }
        }
    }

    /**
     * Generate a permutation matrix that lets us map positions from the precomputed confspace to the new confspace
     */
    public int[] genConfSpaceMapping() {
        // the permutation matrix maps confs in the precomputed flexible to the full confspace
        // Note that I think this works because Positions have equals() check the residue number
        return precomputedPfunc.confSpace.positions.stream()
                .mapToInt(confSpace.positions::indexOf)
                .toArray();
    }

    private class SHARKStarQueue extends PriorityQueue<MultiSequenceSHARKStarNode> {
        private BigDecimal partitionFunctionUpperSum = BigDecimal.ZERO;
        private BigDecimal partitionFunctionLowerSum = BigDecimal.ZERO;
        private final Sequence seq;

        public SHARKStarQueue(Sequence seq) {
            super(new SeqNodeComaparator(seq));
            this.seq = seq;
        }

        public BigDecimal getPartitionFunctionUpperBound() {
            return partitionFunctionUpperSum;
        }

        public BigDecimal getPartitionFunctionLowerBound() {
            return partitionFunctionLowerSum;
        }

        @Override
        public boolean add(MultiSequenceSHARKStarNode node) {
            debugCheck();
            partitionFunctionUpperSum = partitionFunctionUpperSum.add(node.getUpperBound(seq));
            partitionFunctionLowerSum = partitionFunctionLowerSum.add(node.getLowerBound(seq));
            debugCheck();
            return super.add(node);
        }

        @Override
        public MultiSequenceSHARKStarNode poll() {
            MultiSequenceSHARKStarNode node = super.poll();
            debugCheck();
            partitionFunctionUpperSum = partitionFunctionUpperSum.subtract(node.getUpperBound(seq));
            partitionFunctionLowerSum = partitionFunctionLowerSum.subtract(node.getLowerBound(seq));
            debugCheck();
            return node;
        }

        public String toString() {
            return ""+size()+" nodes->"+new MathTools.BigDecimalBounds(partitionFunctionLowerSum, partitionFunctionUpperSum);

        }

        private void debugCheck(boolean force) {
            if(!debug || !force)
                return;
            BigDecimal sumDifference = partitionFunctionUpperSum.subtract(partitionFunctionLowerSum);
            if(sumDifference.compareTo(BigDecimal.ZERO) < 0 && sumDifference.compareTo(BigDecimal.valueOf(1e-5)) > 0)
                System.err.println("Invalid bounds. Lower bound is greater than upper bound.");
            if(partitionFunctionLowerSum.compareTo(BigDecimal.ZERO) < 0)
                System.err.println("Invalid bounds. Lower bound is less than zero.");
            if(!isEmpty() && peek().getLowerBound(seq).compareTo(partitionFunctionLowerSum) > 0)
                System.err.println("The top element is bigger than the entire lower bound sum.");
            assert(sumDifference.compareTo(BigDecimal.ZERO) > 0 || sumDifference.compareTo(BigDecimal.valueOf(1e-5)) <= 0);
            assert (partitionFunctionLowerSum.compareTo(BigDecimal.ZERO) >= 0);
            System.out.println("Queue: bounds "+toString());
            List<MultiSequenceSHARKStarNode> nodes = new ArrayList<>();
            for (int i = 0; i < 10; i++) {
                if(isEmpty())
                    break;
                MultiSequenceSHARKStarNode next = super.poll();
                System.out.println(next.toSeqString(seq));
                nodes.add(next);
            }
            for (MultiSequenceSHARKStarNode node : nodes)
                super.add(node);
        }

        private void debugCheck() {
            debugCheck(false);
        }

    }

    private class SeqNodeComaparator implements Comparator<MultiSequenceSHARKStarNode> {
        private final Sequence seq;
        public SeqNodeComaparator(Sequence seq) {
            this.seq = seq;
        }

        @Override
        public int compare(MultiSequenceSHARKStarNode o1, MultiSequenceSHARKStarNode o2) {
            return -o1.getErrorBound(seq).compareTo(o2.getErrorBound(seq));
        }

        @Override
        public boolean equals(Object obj) {
            return false;
        }

    }

    @Override
    public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
        init(targetEpsilon);
    }

    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        this.status = Status.Estimating;
        if(precomputedPfunc == null)
            precomputeFlexible_expensiveWay();
    }

    /* We are recomputing the energy matrix here because there's some
        stupidly nontrivial parallel array mapping to be done to ensure that
        the extensive energy calculation features we rely on are
        computing the right energy fors the flexible conf space.
     */
    private MultiSequenceSHARKStarBound precomputeFlexible_expensiveWay() {
        SimpleConfSpace flexConfSpace = confSpace.makeFlexibleCopy();
        Sequence unassignedFlex = flexConfSpace.makeUnassignedSequence();
        RCs flexRCs = new RCs(flexConfSpace);

        ConfEnergyCalculator flexMinimizingConfECalc = FlexEmatCalculator_badCode.makeMinimizeConfEcalc(flexConfSpace,
                this.parallelism);
        ConfEnergyCalculator rigidConfECalc = FlexEmatCalculator_badCode.makeRigidConfEcalc(flexMinimizingConfECalc);
        EnergyMatrix flexMinimizingEmat = FlexEmatCalculator_badCode.makeEmat(flexMinimizingConfECalc, "minimizing");
        EnergyMatrix flexRigidEmat = FlexEmatCalculator_badCode.makeEmat(rigidConfECalc, "rigid");
        UpdatingEnergyMatrix flexCorrection = new UpdatingEnergyMatrix(flexConfSpace, flexMinimizingEmat,
                flexMinimizingConfECalc);


        MultiSequenceSHARKStarBound precompFlex = new MultiSequenceSHARKStarBound(
                flexConfSpace, flexRigidEmat, flexMinimizingEmat,
                flexMinimizingConfECalc, flexRCs, this.parallelism);
        precompFlex.initFlex(this.targetEpsilon, this.stabilityThreshold, flexCorrection);
        PartitionFunction flexBound =
                precompFlex.getPartitionFunctionForSequence(unassignedFlex);
        flexBound.compute();
        processPrecomputedFlex(precompFlex);
        return precompFlex;
    }

    public void init(double epsilon, BigDecimal stabilityThreshold) {
        init(epsilon);
        this.stabilityThreshold = stabilityThreshold;
    }

    private void initFlex(double epsilon, BigDecimal stabilityThreshold, UpdatingEnergyMatrix correctionMatrix) {
        this.targetEpsilon = epsilon;
        this.status = Status.Estimating;
        this.stabilityThreshold = stabilityThreshold;
        this.correctionMatrix = correctionMatrix;
    }

    public void setRCs(RCs rcs) {
        fullRCs = rcs;
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



    @Override
    public Status getStatus() {
        return null;
    }

    @Override
    public PartitionFunction.Values getValues() {
        return null;
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
        return numInternalNodesProcessed + numConfsEnergied + numConfsScored + numPartialMinimizations;
    }

    @Override
    public void compute(int maxNumConfs) {
        throw new UnsupportedOperationException("Do not try to run Multisequence SHARK* bounds directly. Call " +
                "makeBoundFor() and use the generated single sequence bound.");
    }

    public void computeForSequence(int maxNumConfs, SingleSequenceSHARKStarBound sequenceBound) {
        debugPrint("Num conformations: " + sequenceBound.numConformations);
        double lastEps = 1;

        int previousConfCount = workDone();

        if (!sequenceBound.nonZeroLower()) {
            runUntilNonZero(sequenceBound);
            sequenceBound.updateBound();
        }

        while (sequenceBound.sequenceEpsilon > targetEpsilon &&
                workDone() - previousConfCount < maxNumConfs
                && isStable(stabilityThreshold, sequenceBound.sequence)) {
            debugPrint("Tightening from epsilon of " + sequenceBound.sequenceEpsilon);
            if (debug) {
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                debugHeap(sequenceBound.fringeNodes);
                printTree(sequenceBound.sequence,rootNode);
            }
            tightenBoundInPhases(sequenceBound);
            debugPrint("Errorbound is now " + sequenceBound.sequenceEpsilon);
            if (lastEps < sequenceBound.sequenceEpsilon && sequenceBound.sequenceEpsilon - lastEps > 0.01
                || sequenceBound.errors()) {
                System.err.println("Error. Bounds got looser.");
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                printTree(sequenceBound.sequence,rootNode);
                System.exit(-1);
            }
            lastEps = sequenceBound.sequenceEpsilon;
        }
        if (!isStable(stabilityThreshold, sequenceBound.sequence))
            sequenceBound.status = Status.Unstable;
        loopTasks.waitForFinish();
        minimizingEcalc.tasks.waitForFinish();
        BigDecimal averageReduction = BigDecimal.ZERO;
        int totalMinimizations = numConfsEnergied + numPartialMinimizations;
        if (totalMinimizations > 0)
            averageReduction = cumulativeZCorrection
                    .divide(new BigDecimal(totalMinimizations), new MathContext(BigDecimal.ROUND_HALF_UP));
        debugPrint(String.format("Average Z reduction per minimization: %12.6e", averageReduction));
    }

    protected void debugPrint(String s) {
        if (debug)
            System.out.println(s);
    }

    protected void profilePrint(String s) {
        if (profileOutput)
            System.out.println(s);
    }

    public void compute() {
        compute(Integer.MAX_VALUE);
    }

    @Override
    public Result makeResult() {
        throw new UnsupportedOperationException("???");
        /*
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
        return result;
        */
    }

    public void setParallelism(Parallelism val) {

        if (val == null) {
            val = Parallelism.makeCpu(1);
        }

        parallelism = val;
        //loopTasks = minimizingEcalc.tasks;
        if (loopTasks == null)
            loopTasks = parallelism.makeTaskExecutor(1000);
        contexts.allocate(parallelism.getParallelism());
    }

    protected boolean shouldMinimize(MultiSequenceSHARKStarNode node, Sequence seq) {
        return node.isLeaf() && !node.isMinimized(seq);
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
        int[] confOfInterest = new int[]{4, 5, 8, 18};
        if (conf.length != confOfInterest.length)
            return;
        boolean match = true;
        for (int i = 0; i < confOfInterest.length; i++) {
            if (conf[i] != confOfInterest[i]) {
                match = false;
                break;
            }
        }
        if (match)
            System.out.println("Matched " + SimpleConfSpace.formatConfRCs(conf));
    }

    // We want to process internal nodes without worrying about the bound too much until we have
    // a nonzero lower bound. We have to have a nonzero lower bound, so we have to have at least
    // one node with a negative conf upper bound.
    private void runUntilNonZero(SingleSequenceSHARKStarBound bound) {
        System.out.println("Running until leaf is found...");
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        double bestConfUpper = Double.POSITIVE_INFINITY;

        List<MultiSequenceSHARKStarNode> newNodes = new ArrayList<>();
        List<MultiSequenceSHARKStarNode> leafNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        if (queue.isEmpty())
            queue.add(rootNode);
        boundLowestBoundConfUnderNode(bound, queue.poll(), newNodes);
        for(MultiSequenceSHARKStarNode newNode : newNodes) {
            if(!newNode.isMinimized(bound.sequence))
                queue.add(newNode);
            else {
                bound.addFinishedNode(newNode);
            }
        }


        newNodes.clear();
        System.out.println("Found a leaf!");
        rootNode.computeEpsilonErrorBounds(bound.sequence);
        //printTree(bound.sequence, rootNode);
        nonZeroLower = true;
    }

    protected void tightenBoundInPhases(SingleSequenceSHARKStarBound bound) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        assert (!queue.isEmpty() || !bound.internalQueue.isEmpty() || !bound.leafQueue.isEmpty());
        System.out.println(String.format("Current overall error bound: %12.10f, spread of [%12.6e, %12.6e]",
                bound.sequenceEpsilon, bound.getValues().calcLowerBound(),
                bound.getValues().calcUpperBound()));
        List<MultiSequenceSHARKStarNode> internalNodes = new ArrayList<>();
        List<MultiSequenceSHARKStarNode> leafNodes = new ArrayList<>();
        List<MultiSequenceSHARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        int numNodes = 0;
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        Stopwatch internalTime = new Stopwatch();
        Stopwatch leafTime = new Stopwatch();
        double leafTimeSum = 0;
        double internalTimeSum = 0;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};
        populateQueues(bound, internalNodes, leafNodes, internalZ, leafZ, ZSums);
        //bound.updateBound();
        debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]", bound.getValues().calcLowerBound(),
                bound.getValues().calcUpperBound()));
        internalZ = ZSums[0];
        leafZ = ZSums[1];
        if(MathTools.isRelativelySame(internalZ, leafZ, PartitionFunction.decimalPrecision, 1e-3)
                && MathTools.isRelativelySame(leafZ, BigDecimal.ZERO, PartitionFunction.decimalPrecision, 1e-3))
            System.out.println("This is a bad time.");
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        if (MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            System.out.println("Processing " + numNodes + " leaf nodes...");
            leafTime.reset();
            leafTime.start();
            for (MultiSequenceSHARKStarNode leafNode : leafNodes) {
                //debugPrint("Processing Node: " + leafNode.toSeqString(bound.sequence));
                processFullConfNode(bound, newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
            }
            loopTasks.waitForFinish();
            leafTime.stop();
            leafTimeAverage = leafTime.getTimeS();
            System.out.println("Processed " + numNodes + " leaves in " + leafTimeAverage + " seconds.");
            if (maxMinimizations < parallelism.numThreads)
                maxMinimizations++;
            bound.internalQueue.addAll(internalNodes);
        } else {
            numNodes = internalNodes.size();
            System.out.println("Processing " + numNodes + " internal nodes...");
            internalTime.reset();
            internalTime.start();
            for (MultiSequenceSHARKStarNode internalNode : internalNodes) {
                if (!MathTools.isGreaterThan(internalNode.getLowerBound(bound.sequence), BigDecimal.ONE) &&
                        MathTools.isGreaterThan(
                                MathTools.bigDivide(internalNode.getUpperBound(bound.sequence), rootNode.getUpperBound(bound.sequence),
                                        PartitionFunction.decimalPrecision),
                                new BigDecimal(1 - targetEpsilon))
                ) {
                    loopTasks.submit(() -> {
                        boundLowestBoundConfUnderNode(bound, internalNode, newNodes);
                        return null;
                    }, (ignored) -> {
                    });
                } else {
                    processPartialConfNode(bound, newNodes, internalNode, internalNode.getConfSearchNode());
                }
                internalNode.markUpdated();
            }
            loopTasks.waitForFinish();
            internalTime.stop();
            internalTimeSum = internalTime.getTimeS();
            internalTimeAverage = internalTimeSum / Math.max(1, internalNodes.size());
            debugPrint("Internal node time :" + internalTimeSum + ", average " + internalTimeAverage);
            numInternalNodesProcessed += internalNodes.size();
            bound.leafQueue.addAll(leafNodes);
        }
        loopCleanup(bound, newNodes, loopWatch, numNodes);
    }

    protected void debugHeap(SHARKStarQueue queue) {
        queue.debugCheck(true);
    }


    boolean isStable(BigDecimal stabilityThreshold, Sequence seq) {
        return numConfsEnergied <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(rootNode.getUpperBound(seq), stabilityThreshold);
    }


    protected void populateQueues(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> internalNodes, List<MultiSequenceSHARKStarNode> leafNodes, BigDecimal internalZ,
                                  BigDecimal leafZ, BigDecimal[] ZSums) {
        List<MultiSequenceSHARKStarNode> leftoverLeaves = new ArrayList<>();
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        int maxNodes = 1000;
        //int maxNodes = 1;
        if (leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int) Math.floor(0.1 * leafTimeAverage / internalTimeAverage));
        while (!queue.isEmpty() && (bound.internalQueue.size() < maxNodes || bound.leafQueue.size() < maxMinimizations)) {
            MultiSequenceSHARKStarNode curNode = queue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(fullRCs.getNumPos());

            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore =  -node.getPartialConfLowerBound();
            double confCorrection = Math.min(correctgscore, node.getPartialConfUpperBound()) + hscore;
            if (!curNode.isMinimized(bound.sequence) && curNode.getConfLowerBound(bound.sequence) < confCorrection
                    && curNode.getConfLowerBound(bound.sequence) - confCorrection > 1e-5) {
                if (confCorrection < curNode.getConfLowerBound(bound.sequence)) {
                    System.out.println("huh!?");
                }
                System.out.println("Correction from " + correctionMatrix.sourceECalc + ":" + node.getPartialConfLowerBound() + "->" + correctgscore);
                recordCorrection(curNode.getConfLowerBound(bound.sequence), correctgscore - node.getPartialConfLowerBound());

                node.setPartialConfLowerAndUpper(correctgscore, node.getPartialConfUpperBound());
                if (confCorrection > node.getPartialConfUpperBound()) {
                    System.out.println("Overcorrected" + SimpleConfSpace.formatConfRCs(node.assignments) + ": " + confCorrection + " > " + node.getPartialConfUpperBound());
                    node.setPartialConfLowerAndUpper(node.getPartialConfUpperBound(), node.getPartialConfUpperBound());
                    confCorrection = node.getPartialConfUpperBound() + hscore;
                }
                node.setBoundsFromConfLowerAndUpper(confCorrection, curNode.getConfUpperBound(bound.sequence));
                curNode.markUpdated();
                leftoverLeaves.add(curNode);
                continue;
            }


            if (node.getLevel() < fullRCs.getNumPos()) {
                bound.internalQueue.add(curNode);
            } else if (shouldMinimize(curNode, bound.sequence) && !correctedNode(leftoverLeaves, curNode, node, bound.sequence)) {
                bound.leafQueue.add(curNode);
            }

        }

        ZSums[0] = fillListFromQueue(internalNodes, bound.internalQueue, maxNodes, bound.sequence);
        ZSums[1] = fillListFromQueue(leafNodes, bound.leafQueue, maxMinimizations, bound.sequence);
        queue.addAll(leftoverLeaves);
    }

    private BigDecimal fillListFromQueue(List<MultiSequenceSHARKStarNode> list, Queue<MultiSequenceSHARKStarNode> queue, int max, Sequence seq) {
        BigDecimal sum = BigDecimal.ZERO;
        List<MultiSequenceSHARKStarNode> leftovers = new ArrayList<>();
        while (!queue.isEmpty() && list.size() < max) {
            MultiSequenceSHARKStarNode curNode = queue.poll();
            if (correctedNode(leftovers, curNode, curNode.getConfSearchNode(), seq)) {
                continue;
            }
            BigDecimal diff = curNode.getUpperBound(seq).subtract(curNode.getLowerBound(seq));
            sum = sum.add(diff);
            list.add(curNode);
        }
        queue.addAll(leftovers);
        return sum;
    }

    protected void loopCleanup(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> newNodes, Stopwatch loopWatch, int numNodes) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        for (MultiSequenceSHARKStarNode node : newNodes) {
            if (node != null){
                if(node.isMinimized(bound.sequence))
                    bound.addFinishedNode(node);
                else
                    queue.add(node);
            }
        }
        bound.updateBound();
        loopWatch.stop();
        double loopTime = loopWatch.getTimeS();
        profilePrint("Processed " + numNodes + " this loop, spawning " + newNodes.size() + " in " + loopTime + ", " + stopwatch.getTime() + " so far");
        loopWatch.reset();
        loopWatch.start();
        processPreminimization(bound, minimizingEcalc);
        profilePrint("Preminimization time : " + loopWatch.getTime(2));
        double curEpsilon = bound.sequenceEpsilon;
        //printTree(bound.sequence, rootNode);
        //rootNode.updateConfBounds(new ConfIndex(RCs.getNumPos()), RCs, gscorer, hscorer);
        loopWatch.stop();
        cleanupTime = loopWatch.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
        System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]", bound.getValues().calcLowerBound(),
                bound.getValues().calcUpperBound()));
    }

    protected boolean correctedNode(List<MultiSequenceSHARKStarNode> newNodes, MultiSequenceSHARKStarNode curNode, Node node, Sequence seq) {
        assert (curNode != null && node != null);
        double confCorrection = correctionMatrix.confE(node.assignments);
        if (node.getPartialConfLowerBound() < confCorrection) {
            double oldg = node.getPartialConfLowerBound();
            node.setBoundsFromConfLowerAndUpper(confCorrection, node.getPartialConfUpperBound());
            recordCorrection(oldg, confCorrection - oldg);
            //node.setBoundsFromConfLowerAndUpper(curNode.getConfLowerBound(seq) - oldg + confCorrection, curNode.getConfUpperBound(seq));
            curNode.setBoundsFromConfLowerAndUpper(curNode.getConfLowerBound(seq) - oldg + confCorrection, curNode.getConfUpperBound(seq), seq);
            curNode.markUpdated();
            newNodes.add(curNode);
            return true;
        }
        return false;
    }

    private MultiSequenceSHARKStarNode drillDown(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> newNodes,
                                                 MultiSequenceSHARKStarNode curNode, Node node) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            node.index(context.index);
            // which pos to expand next?
            int nextPos = order.getNextPos(context.index, RCs);
            assert (!context.index.isDefined(nextPos));
            assert (context.index.isUndefined(nextPos));

            // score child nodes with tasks (possibly in parallel)
            List<MultiSequenceSHARKStarNode> children = new ArrayList<>();
            double bestChildLower = Double.POSITIVE_INFINITY;
            MultiSequenceSHARKStarNode bestChild = null;
            for (int nextRc : RCs.get(nextPos)) {

                if (hasPrunedPair(context.index, nextPos, nextRc)) {
                    continue;
                }

                /** We don't currently prune. To do so, we'd need to preserve some code we don't use. */
                // if this child was pruned dynamically, then don't score it
                /*
                if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                    continue;
                }
                */
                Stopwatch partialTime = new Stopwatch().start();
                Node child = node.assign(nextPos, nextRc);
                double confLowerBound = Double.POSITIVE_INFINITY;
                double confUpperBound = Double.NEGATIVE_INFINITY;

                if(confMatch(child.assignments, debugConf))
                    System.out.println("Gotcha-prequel-drill.");

                // score the child node differentially against the parent node
                if (child.getLevel() < RCs.getNumPos()) {
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double hdiff = context.lowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double maxhdiff = context.upperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    //Correct for incorrect gscore.
                    //rigiddiff = rigiddiff - node.partialConfLowerbound + node.partialConfUpperBound;
                    child.setPartialConfLowerAndUpper(diff, rigiddiff);

                    confLowerBound = child.getPartialConfLowerBound() + hdiff;
                    confUpperBound = rigiddiff + maxhdiff;
                    child.computeNumConformations(RCs);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound);
                    progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.sequenceEpsilon);
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    //confRigid = confRigid - node.partialConfLowerbound + node.partialConfUpperBound;

                    child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                    double confCorrection = correctionMatrix.confE(child.assignments);
                    double lowerbound = minimizingEmat.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    checkBounds(confCorrection, confRigid);
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    confLowerBound = lowerbound;
                    child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                    confUpperBound = confRigid;
                    numConfsScored++;
                    progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.sequenceEpsilon);
                }
                partialTime.stop();
                loopPartialTime += partialTime.getTimeS();


                MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild = curNode.makeChild(child, bound.sequence, confLowerBound, confUpperBound);
                if (Double.isNaN(child.getPartialConfUpperBound()))
                    System.out.println("Huh!?");
                MultiSequenceSHARKStarNodeChild.markUpdated();
                if (confLowerBound < bestChildLower) {
                    bestChild = MultiSequenceSHARKStarNodeChild;
                }
                // collect the possible children
                if (MultiSequenceSHARKStarNodeChild.getConfLowerBound(bound.sequence) < 0) {
                    children.add(MultiSequenceSHARKStarNodeChild);
                }
                newNodes.add(MultiSequenceSHARKStarNodeChild);

            }
            return bestChild;
        }
    }

    protected void boundLowestBoundConfUnderNode(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarNode startNode,
                                                 List<MultiSequenceSHARKStarNode> generatedNodes) {
        System.out.println("Bounding "+startNode.getConfSearchNode());
        Comparator<MultiSequenceSHARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getConfLowerBound(bound.sequence));
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        PriorityQueue<MultiSequenceSHARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
        drillQueue.add(startNode);

        List<MultiSequenceSHARKStarNode> newNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        while (!drillQueue.isEmpty()) {
            numNodes++;
            MultiSequenceSHARKStarNode curNode = drillQueue.poll();
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(RCs.getNumPos());
            node.index(index);

            if (node.getLevel() < RCs.getNumPos()) {
                MultiSequenceSHARKStarNode nextNode = drillDown(bound, newNodes, curNode, node);
                newNodes.remove(nextNode);
                drillQueue.add(nextNode);
            } else {
                newNodes.add(curNode);
            }

            //debugHeap(drillQueue, true);
            if (leafLoop.getTimeS() > 1) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",
                        numNodes,
                        overallLoop.getTime(2),
                        rootNode.getLowerBound(bound.sequence),
                        rootNode.getUpperBound(bound.sequence)));
            }
        }
        generatedNodes.addAll(newNodes);

    }

    private class PartialConfNodeResult {
        Node resultNode = null;
        double upperBound = Double.NaN;
        double lowerBound = Double.NaN;

        public boolean isValid() {
            return resultNode != null && !Double.isNaN(upperBound) && !Double.isNaN(lowerBound);
        }
    }

    protected void processPartialConfNode(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> newNodes,
                                          MultiSequenceSHARKStarNode curNode, Node node) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        //debugPrint("Processing "+curNode.toSeqString(bound.sequence));
        // which pos to expand next?
        node.index(confIndex);
        int nextPos = order.getNextPos(confIndex, RCs);
        assert (!confIndex.isDefined(nextPos));
        assert (confIndex.isUndefined(nextPos));
        if(confMatch(curNode.getConfSearchNode().assignments, debugConf)) {
            System.out.println("Gotcha.");
            System.out.println("Energy of partial conf: "+rigidEmat.getInternalEnergy(curNode.toTuple()));
        }

        // score child nodes with tasks (possibly in parallel)
        List<MultiSequenceSHARKStarNode> children = new ArrayList<>();
        for (int nextRc : bound.seqRCs.get(nextPos)) {

            if (hasPrunedPair(confIndex, nextPos, nextRc)) {
                continue;
            }

            /** We don't currently prune. To do so, we'd need to preserve some code we don't use. */
            // if this child was pruned dynamically, then don't score it
            /*
            if (pruner != null && pruner.isPruned(node, nextPos, nextRc)) {
                continue;
            }
            */

            loopTasks.submit(() -> {

                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    Stopwatch partialTime = new Stopwatch().start();
                    ScoreContext context = checkout.get();
                    node.index(context.index);
                    Node child = node.assign(nextPos, nextRc);
                    PartialConfNodeResult result = new PartialConfNodeResult();
                    result.resultNode = child;
                    if(confMatch(child.assignments, debugConf))
                        System.out.println("Gotcha-prequel.");

                    // score the child node differentially against the parent node
                    if (child.getLevel() < RCs.getNumPos()) {
                        double wtf;
                        if(confMatch(child.assignments, debugConf)) {
                            wtf = rigidEmat.getInternalEnergy(curNode.toTuple());
                            double wtff = curNode.getConfSearchNode().getRigidGScore();
                            double diffwtf = wtf-wtff;
                            System.out.println("Difftf? "+diffwtf);
                        }
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double diff = confCorrection;
                        double rigiddiff = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double hdiff = context.lowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double maxhdiff = context.upperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double checkNum = 0;
                        if(confMatch(child.assignments, debugConf)) {
                            checkNum = rigidEmat.getInternalEnergy(new RCTuple(child.makeConf(fullRCs.getNumPos())));
                            double diffNum = rigiddiff - checkNum;
                            ConfIndex copy = context.index.assign(nextPos, nextRc);
                            double checkNum2 = context.partialConfUpperBoundScorer.calc(copy, RCs);
                            System.out.println("The weirdness is "+diffNum);
                        }
                        //Correct for incorrect gscore.
                        //rigiddiff = rigiddiff - node.partialConfLowerbound + node.partialConfUpperBound;
                        child.setPartialConfLowerAndUpper(diff, rigiddiff);

                        double confLowerBound = child.getPartialConfLowerBound() + hdiff;
                        double confUpperbound = rigiddiff + maxhdiff;
                        child.computeNumConformations(RCs);
                        double lowerbound = minimizingEmat.confE(child.assignments);
                        if (diff < confCorrection) {
                            recordCorrection(confLowerBound, confCorrection - diff);
                            confLowerBound = confCorrection + hdiff;
                        }
                        child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);
                        progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.sequenceEpsilon);
                        result.lowerBound = confLowerBound;
                        result.upperBound = confUpperbound;
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        //confRigid = confRigid - node.partialConfLowerbound + node.partialConfUpperBound;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = minimizingEmat.confE(child.assignments);

                        if (lowerbound < confCorrection) {
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection, confRigid);
                        child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                        child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                        numConfsScored++;
                        progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.sequenceEpsilon);
                        result.lowerBound = confCorrection;
                        result.upperBound = confRigid;
                    }
                    partialTime.stop();
                    loopPartialTime += partialTime.getTimeS();


                    return result;
                }

            }, (PartialConfNodeResult result) -> {
                assert (result.isValid());
                MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild = curNode.makeChild(result.resultNode,
                        bound.sequence, result.lowerBound, result.upperBound);
                MultiSequenceSHARKStarNodeChild.setBoundsFromConfLowerAndUpper(result.lowerBound,
                        result.upperBound, bound.sequence);
                //System.out.println("Created new child "+MultiSequenceSHARKStarNodeChild.toSeqString(bound.sequence));
                // collect the possible children
                if (MultiSequenceSHARKStarNodeChild.getConfLowerBound(bound.sequence) < 0) {
                    children.add(MultiSequenceSHARKStarNodeChild);

                }
                if (!MultiSequenceSHARKStarNodeChild.isMinimized(bound.sequence)) {
                    newNodes.add(MultiSequenceSHARKStarNodeChild);
                } else {
                    MultiSequenceSHARKStarNodeChild.computeEpsilonErrorBounds(bound.sequence);
                    bound.addFinishedNode(MultiSequenceSHARKStarNodeChild);
                }

                //curNode.updateSubtreeBounds(bound.sequence);
                //printTree(bound.sequence, curNode);
                curNode.markUpdated();
            });
        }
    }

    private boolean confMatch(int[] assignments, int[] ints) {
        if(assignments.length != ints.length)
            return false;
        for(int i = 0; i < assignments.length; i++) {
           if(assignments[i] != ints[i])
               return false;
        }
        return true;
    }


    protected void processFullConfNode(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> newNodes,
                                       MultiSequenceSHARKStarNode curNode, Node node) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        double confCorrection = correctionMatrix.confE(node.assignments);
        if (curNode.getConfLowerBound(bound.sequence) < confCorrection || node.getPartialConfLowerBound() < confCorrection) {
            double oldg = node.getPartialConfLowerBound();
            node.setPartialConfLowerAndUpper(confCorrection, confCorrection);
            recordCorrection(oldg, confCorrection - oldg);
            curNode.setBoundsFromConfLowerAndUpper(confCorrection, curNode.getConfUpperBound(bound.sequence), bound.sequence);
            node.setBoundsFromConfLowerAndUpper(confCorrection, curNode.getConfUpperBound(bound.sequence));
            curNode.markUpdated();
            newNodes.add(curNode);
            System.out.println("Correcting "+curNode.toSeqString(bound.sequence));
            return;
        }
        loopTasks.submit(() -> {
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        ScoreContext context = checkout.get();
                        node.index(context.index);

                        System.out.println("Minmized "+curNode.toSeqString(bound.sequence));
                        ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, curNode.getConfLowerBound(bound.sequence));
                        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                        Stopwatch correctionTimer = new Stopwatch().start();
                        computeEnergyCorrection(analysis, conf, context.ecalc, bound.sequenceEpsilon);

                        double energy = analysis.epmol.energy;
                        double newConfUpper = energy;
                        double newConfLower = energy;
                        // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
                        double oldConfUpper = curNode.getConfUpperBound(bound.sequence);
                        double oldConfLower = curNode.getConfLowerBound(bound.sequence);
                        if (newConfUpper > oldConfUpper) {
                            System.err.println("Upper bounds got worse after minimization:" + newConfUpper
                                    + " > " + (oldConfUpper) + ". Rejecting minimized energy.");
                            System.err.println("Node info: " + node);

                            newConfUpper = oldConfUpper;
                            newConfLower = oldConfUpper;
                        }
                        curNode.setBoundsFromConfLowerAndUpper(newConfLower, newConfUpper, bound.sequence);
                        double oldgscore = node.getPartialConfLowerBound();
                        node.setPartialConfLowerAndUpper(newConfLower, newConfUpper);
                        String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (curNode.getConfLowerBound(bound.sequence)) + "," + (curNode.getConfUpperBound(bound.sequence)) + "]";
                        debugPrint(out);
                        curNode.markUpdated();
                        synchronized (this) {
                            numConfsEnergied++;
                            minList.set(conf.getAssignments().length - 1, minList.get(conf.getAssignments().length - 1) + 1);
                            recordReduction(oldConfLower, oldConfUpper, energy);
                            printMinimizationOutput(node, newConfLower, oldgscore, bound);
                            bound.addFinishedNode(curNode);
                        }


                    }
                    return null;
                },
                // Dummy function. We're not doing anything here.
                (Node child) -> {
                    progress.reportLeafNode(node.getPartialConfLowerBound(), queue.size(), bound.sequenceEpsilon);
                    if (!curNode.isMinimized(bound.sequence))
                        newNodes.add(curNode);
                });
    }

    private void printMinimizationOutput(Node node, double newConfLower, double oldgscore, SingleSequenceSHARKStarBound seqBound) {
        if (printMinimizedConfs) {
            System.out.println("[" + SimpleConfSpace.formatConfRCs(node.assignments) + "]"
                    + String.format("conf:%4d, score:%12.6f, lower:%12.6f, corrected:%12.6f energy:%12.6f"
                            + ", bounds:[%12e, %12e], delta:%12.6f, time:%10s",
                    numConfsEnergied, oldgscore, minimizingEmat.confE(node.assignments),
                    correctionMatrix.confE(node.assignments), newConfLower,
                    rootNode.getLowerBound(seqBound.sequence), rootNode.getUpperBound(seqBound.sequence),
                    seqBound.sequenceEpsilon, stopwatch.getTime(2)));

        }
    }


    private void checkBounds(double lower, double upper) {
        if (upper < lower && upper - lower > 1e-5 && upper < 10)
            debugPrint("Bounds incorrect.");
    }

    private void computeEnergyCorrection(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                         ConfEnergyCalculator ecalc, double epsilonBound) {
        if (conf.getAssignments().length < 3)
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
        for (int pos = 0; pos < diff.getNumPos(); pos++) {
            for (int rc = 0; rc < diff.getNumConfAtPos(pos); rc++) {
                for (int pos2 = 0; pos2 < diff.getNumPos(); pos2++) {
                    for (int rc2 = 0; rc2 < diff.getNumConfAtPos(pos2); rc2++) {
                        if (pos >= pos2)
                            continue;
                        double sum = 0;
                        sum += diff.getOneBody(pos, rc);
                        sum += diff.getPairwise(pos, rc, pos2, rc2);
                        sum += diff.getOneBody(pos2, rc2);
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
        for (int i = 0; i < sortedPairwiseTerms2.size(); i++) {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if (pairDiff < minDifference && maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            int localMinimizations = 0;
            for (int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double tupleBounds = rigidEmat.getInternalEnergy(tuple) - minimizingEmat.getInternalEnergy(tuple);
                if (tupleBounds < triplethreshhold)
                    continue;
                minList.set(tuple.size() - 1, minList.get(tuple.size() - 1) + 1);
                computeDifference(tuple, minimizingEcalc);
                localMinimizations++;
            }
            numPartialMinimizations += localMinimizations;
            progress.reportPartialMinimization(localMinimizations, epsilonBound);
        }
        correctionTime.stop();
        ecalc.tasks.waitForFinish();
    }


    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if (correctedTuples.contains(tuple.stringListing()))
            return;
        correctedTuples.add(tuple.stringListing());
        if (correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
            double tripleEnergy = minimizedTuple.energy;

            double lowerbound = minimizingEmat.getInternalEnergy(tuple);
            if (tripleEnergy - lowerbound > 0) {
                double correction = tripleEnergy - lowerbound;
                correctionMatrix.setHigherOrder(tuple, correction);
            } else
                System.err.println("Negative correction for " + tuple.stringListing());
        });
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for (int pos : positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }

    private void processPreminimization(SingleSequenceSHARKStarBound bound, ConfEnergyCalculator ecalc) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        int maxMinimizations = 1;//parallelism.numThreads;
        List<MultiSequenceSHARKStarNode> topConfs = getTopConfs(queue, maxMinimizations);
        // Need at least two confs to do any partial preminimization
        if (topConfs.size() < 2) {
            queue.addAll(topConfs);
            return;
        }
        RCTuple lowestBoundTuple = topConfs.get(0).toTuple();
        RCTuple overlap = findLargestOverlap(lowestBoundTuple, topConfs, 3);
        //Only continue if we have something to minimize
        for (MultiSequenceSHARKStarNode conf : topConfs) {
            RCTuple confTuple = conf.toTuple();
            if (minimizingEmat.getInternalEnergy(confTuple) == rigidEmat.getInternalEnergy(confTuple))
                continue;
            numPartialMinimizations++;
            minList.set(confTuple.size() - 1, minList.get(confTuple.size() - 1) + 1);
            if (confTuple.size() > 2 && confTuple.size() < RCs.getNumPos()) {
                minimizingEcalc.tasks.submit(() -> {
                    computeTupleCorrection(minimizingEcalc, conf.toTuple(), bound.sequenceEpsilon);
                    return null;
                }, (econf) -> {
                });
            }
        }
        //minimizingEcalc.tasks.waitForFinish();
        if (overlap.size() > 3 && !correctionMatrix.hasHigherOrderTermFor(overlap)
                && minimizingEmat.getInternalEnergy(overlap) != rigidEmat.getInternalEnergy(overlap)) {
            minimizingEcalc.tasks.submit(() -> {
                computeTupleCorrection(ecalc, overlap, bound.sequenceEpsilon);
                return null;
            }, (econf) -> {
            });
        }
        queue.addAll(topConfs);
    }

    private void computeTupleCorrection(ConfEnergyCalculator ecalc, RCTuple overlap, double epsilonBound) {
        if (correctionMatrix.hasHigherOrderTermFor(overlap))
            return;
        double pairwiseLower = minimizingEmat.getInternalEnergy(overlap);
        double partiallyMinimizedLower = ecalc.calcEnergy(overlap).energy;
        progress.reportPartialMinimization(1, epsilonBound);
        if (partiallyMinimizedLower > pairwiseLower)
            synchronized (correctionMatrix) {
                correctionMatrix.setHigherOrder(overlap, partiallyMinimizedLower - pairwiseLower);
            }
        progress.reportPartialMinimization(1, epsilonBound);
    }

    private List<MultiSequenceSHARKStarNode> getTopConfs(PriorityQueue<MultiSequenceSHARKStarNode> queue, int numConfs) {
        List<MultiSequenceSHARKStarNode> topConfs = new ArrayList<>();
        while (topConfs.size() < numConfs && !queue.isEmpty()) {
            MultiSequenceSHARKStarNode nextLowestConf = queue.poll();
            topConfs.add(nextLowestConf);
        }
        return topConfs;
    }


    private RCTuple findLargestOverlap(RCTuple conf, List<MultiSequenceSHARKStarNode> otherConfs, int minResidues) {
        RCTuple overlap = conf;
        for (MultiSequenceSHARKStarNode other : otherConfs) {
            overlap = overlap.intersect(other.toTuple());
            if (overlap.size() < minResidues)
                break;
        }
        return overlap;

    }

    private boolean hasPrunedPair(ConfIndex confIndex, int nextPos, int nextRc) {

        // do we even have pruned pairs?
        PruningMatrix pmat = fullRCs.getPruneMat();
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

    public void setCorrections(UpdatingEnergyMatrix cachedCorrections) {
        correctionMatrix = cachedCorrections;
    }


    public static class Values extends PartitionFunction.Values {

        public Values() {
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


    protected static class ScoreContext {
        public ConfIndex index;
        public AStarScorer partialConfLowerBoundScorer;
        public AStarScorer lowerBoundScorer;
        public AStarScorer upperBoundScorer;
        public AStarScorer partialConfUpperBoundScorer;
        public ConfEnergyCalculator ecalc;
    }

    private class SHARKStarNodeScorer implements AStarScorer {

        private EnergyMatrix emat;
        private MathTools.Optimizer opt = MathTools.Optimizer.Minimize;


        public SHARKStarNodeScorer(EnergyMatrix emat) {
            this(emat, false);
        }
        public SHARKStarNodeScorer(EnergyMatrix emat, boolean negated) {
            this.emat = emat;
            if(negated)
                opt = MathTools.Optimizer.Maximize;
        }

        public SHARKStarNodeScorer(EnergyMatrix emats, RCs rcs) {
           this(emats);
        }

        public SHARKStarNodeScorer(EnergyMatrix emats, RCs rcs, boolean negated) {
            this(emats, negated);
        }

        @Override
        public AStarScorer make() {
            return new SHARKStarNodeScorer(emat);
        }


        public double calc(ConfIndex confIndex, Sequence seq, SimpleConfSpace confSpace) {
            return calc(confIndex, seq.makeRCs(confSpace));
        }

        /* Assumes: that the rcs contain only the sequence in question. In this case, we need only
         *  sum over all unassigned positions. Returns a lower bound on the ensemble energy.
         *  Note: I currently exponentiate and log for compatibilty. This could be optimized.*/
        @Override
        public double calc(ConfIndex confIndex, edu.duke.cs.osprey.astar.conf.RCs rcs) {
            BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
            BigDecimal pfuncBound = BigDecimal.ONE;
            for (int undefinedPosIndex1 = 0; undefinedPosIndex1 < confIndex.numUndefined; undefinedPosIndex1++) {
                int undefinedPos1 = confIndex.undefinedPos[undefinedPosIndex1];
                BigDecimal residueSum = BigDecimal.ZERO;
                for (int rot1 : rcs.get(undefinedPos1)) {
                    double rotEnergy = emat.getEnergy(undefinedPos1, rot1);
                    for (int definedPosIndex = 0; definedPosIndex < confIndex.numDefined; definedPosIndex ++) {
                        int definedPos = confIndex.definedPos[definedPosIndex];
                        int definedRC = confIndex.definedRCs[definedPosIndex];
                        rotEnergy += emat.getEnergy(undefinedPos1, rot1, definedPos, definedRC);
                    }
                    for (int undefinedPosIndex2 = 0; undefinedPosIndex2 < confIndex.numUndefined; undefinedPosIndex2++) {
                        int undefinedPos2 = confIndex.undefinedPos[undefinedPosIndex2];
                        if (undefinedPos2 >= undefinedPos1)
                            continue;
                        double bestPair = Double.MAX_VALUE;
                        for (int rot2 : rcs.get(undefinedPos2)) {
                            bestPair = opt.opt(bestPair, emat.getEnergy(undefinedPos1, rot1, undefinedPos2, rot2));
                        }
                        rotEnergy+= bestPair;
                    }
                    int[] conf = confIndex.makeConf();
                    if(confMatch(conf, debugConf)){
                        System.out.println("Gotcha-calc");
                        conf[undefinedPos1] = rot1;
                        RCTuple testTuple = new RCTuple(conf);
                        double energyCheck = emat.getInternalEnergy(testTuple);
                        System.out.println("Energy of " + testTuple + ": " + energyCheck);
                        System.out.println("Rot energy of "+ testTuple + ": " + rotEnergy);
                    }
                    residueSum = residueSum.add(bcalc.calc(rotEnergy));
                }
                int[] conf = confIndex.makeConf();
                if(confMatch(conf, debugConf)){
                    System.out.println("Gotcha-calc2");
                    System.out.println("End residue sum: "+residueSum);
                }
                pfuncBound = pfuncBound.multiply(residueSum, PartitionFunction.decimalPrecision);
            }
            return bcalc.freeEnergy(pfuncBound);
        }

    }

    /**
     * Thin wrapper class to play nice with BBK* and MSK*
     */
    public class SingleSequenceSHARKStarBound implements PartitionFunction {
        public final Sequence sequence;
        private MultiSequenceSHARKStarBound multisequenceBound;
        private Status status;
        private MultiSequenceSHARKStarBound.Values values;
        private int numConfsEvaluated = 0;
        public final BigInteger numConformations;
        private SHARKStarQueue fringeNodes;
        public SHARKStarQueue internalQueue;
        public SHARKStarQueue leafQueue;
        private double sequenceEpsilon = 1;
        private BigDecimal finishedNodeZ = BigDecimal.ZERO;
        public final RCs seqRCs;

        //debug variable
        private Set<MultiSequenceSHARKStarNode> finishedNodes = new HashSet<>();
        private boolean errors;

        public SingleSequenceSHARKStarBound(Sequence seq, MultiSequenceSHARKStarBound sharkStarBound) {
            this.sequence = seq;
            this.multisequenceBound = sharkStarBound;
            this.seqRCs = seq.makeRCs(sharkStarBound.confSpace);
            this.numConformations = seqRCs.getNumConformations();
            this.fringeNodes = new SHARKStarQueue(seq);
            this.internalQueue = new SHARKStarQueue(seq);
            this.leafQueue = new SHARKStarQueue(seq);
        }

        public void addFinishedNode(MultiSequenceSHARKStarNode node) {
            finishedNodeZ = finishedNodeZ.add(node.getUpperBound(sequence));
            System.out.println("Adding "+node.toSeqString(sequence)+" to finished set");
            if(finishedNodes.contains(node))
                System.err.println("Dupe node addition.");
            finishedNodes.add(node);
        }

        @Override
        public void setReportProgress(boolean val) {
            multisequenceBound.setReportProgress(val);
        }

        @Override
        public void setConfListener(ConfListener val) {
            multisequenceBound.setConfListener(val);
        }

        @Override
        public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
            init(confSearch, null, numConfsBeforePruning, targetEpsilon);
        }

        @Override
        public void init(ConfSearch upperBoundConfs, ConfSearch lowerBoundConfs, BigInteger numConfsBeforePruning, double targetEpsilon) {
            values = new MultiSequenceSHARKStarBound.Values();
            status = Status.Estimating;
        }


        @Override
        public void setStabilityThreshold(BigDecimal stabilityThreshold) {
            multisequenceBound.setStabilityThreshold(stabilityThreshold);
        }

        @Override
        public Status getStatus() {
            return this.status;
        }

        @Override
        public Values getValues() {
            return this.values;
        }

        @Override
        public int getParallelism() {
            return multisequenceBound.getParallelism();
        }

        @Override
        public int getNumConfsEvaluated() {
            return numConfsEvaluated;
        }

        @Override
        public void compute(int maxNumConfs) {
            multisequenceBound.computeForSequence(maxNumConfs, this);
            updateBound();
            if (sequenceEpsilon < targetEpsilon) {
                status = Status.Estimated;
                if (values.qstar.compareTo(BigDecimal.ZERO) == 0) {
                    status = Status.Unstable;
                }
            }
        }

        @Override
        public void compute() {
            compute(Integer.MAX_VALUE);
        }

        @Override
        public Result makeResult() {
            // Calculate the upper bound z reductions from conf lower bounds, since we don't explicitly record these
            lowerReduction_ConfUpperBound = rootNode.getLowerBound(sequence)
                    .subtract(startLowerBound).subtract(lowerReduction_FullMin);
            // Calculate the lower bound z reductions from conf upper bounds, since we don't explicitly record these
            upperReduction_ConfLowerBound = startUpperBound.subtract(rootNode.getUpperBound(sequence))
                    .subtract(upperReduction_FullMin).subtract(upperReduction_PartialMin);

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

        public void updateBound() {
            rootNode.computeEpsilonErrorBounds(sequence);
            BigDecimal lastUpper = values.pstar;
            BigDecimal lastLower = values.qstar;
            BigDecimal upperBound = fringeNodes.getPartitionFunctionUpperBound()
                    .add(internalQueue.getPartitionFunctionUpperBound())
                    .add(leafQueue.getPartitionFunctionUpperBound())
                    .add(finishedNodeZ);
            BigDecimal lowerBound = fringeNodes.getPartitionFunctionLowerBound()
                    .add(internalQueue.getPartitionFunctionLowerBound())
                    .add(leafQueue.getPartitionFunctionLowerBound())
                    .add(finishedNodeZ);
            if(MathTools.isLessThan(lowerBound, lastLower) && ! MathTools.isRelativelySame(lowerBound, lastLower, PartitionFunction.decimalPrecision, 10)) {
                System.err.println("Bounds getting looser. Lower bound is getting lower...");
                errors = true;
            }
            if(MathTools.isGreaterThan(upperBound, lastUpper) && ! MathTools.isRelativelySame(upperBound, lastUpper, PartitionFunction.decimalPrecision, 10)) {
                System.err.println("Bounds getting looser. Upper bound is getting bigger...");
                errors = true;
            }
            values.pstar = upperBound;
            values.qstar = lowerBound;
            values.qprime = upperBound;
            if (upperBound.subtract(lowerBound).compareTo(BigDecimal.ONE) < 1) {
                sequenceEpsilon = 0;
            } else {
                sequenceEpsilon = upperBound.subtract(lowerBound)
                        .divide(upperBound, RoundingMode.HALF_UP).doubleValue();
            }
        }
        List<SeqSpace.ResType> getRTs(SimpleConfSpace.Position confPos, SeqAStarNode.Assignments assignments) {

            // TODO: pre-compute this somehow?
            SeqSpace seqSpace = sequence.seqSpace;

            // map the conf pos to a sequence pos
            SeqSpace.Position seqPos = seqSpace.getPosition(confPos.resNum);
            if (seqPos != null) {

                Integer assignedRT = assignments.getAssignment(seqPos.index);
                if (assignedRT != null) {
                    // use just the assigned res type
                    return Collections.singletonList(seqPos.resTypes.get(assignedRT));
                } else {
                    // use all the res types at the pos
                    return seqPos.resTypes;
                }

            } else {

                // immutable position, use all the res types (should just be one)
                assert (confPos.resTypes.size() == 1);

                // use the null value to signal there's no res type here
                return Collections.singletonList(null);
            }
        }

        List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, SeqSpace.ResType rt, SHARKStar.State state) {
            // TODO: pre-compute this somehow?
            if (rt != null) {
                // mutable pos, grab the RCs that match the RT
                return pos.resConfs.stream()
                        .filter(rc -> rc.template.name.equals(rt.name))
                        .collect(Collectors.toList());
            } else {
                // immutable pos, use all the RCs
                return pos.resConfs;
            }
        }

        public boolean nonZeroLower() {
            return this.fringeNodes.getPartitionFunctionLowerBound().compareTo(BigDecimal.ZERO) > 0;
        }

        public boolean errors() {
            return errors;
        }
    }

    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }

    private static class FlexEmatCalculator_badCode {

        public static EnergyMatrix makeEmat(ConfEnergyCalculator confECalc, String name) {
            System.out.println("Making energy matrix for "+confECalc);
            EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confECalc)
                    .setCacheFile(new File("flex"+"."+name+".emat"))
                    .build()
                    .calcEnergyMatrix();
            return emat;
        }

        public static ConfEnergyCalculator makeRigidConfEcalc(ConfEnergyCalculator sourceConfEcalc){
            EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(sourceConfEcalc.ecalc)
                    .setIsMinimizing(false)
                    .build();
            ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(sourceConfEcalc, ecalcRigid);
            return confEcalcRigid;
        }

        public static ConfEnergyCalculator makeMinimizeConfEcalc(SimpleConfSpace confSpace, Parallelism parallelism) {
            ForcefieldParams ffparams = new ForcefieldParams();

            // how should we compute energies of molecules?
            EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                    .setParallelism(parallelism)
                    .build();
            // how should we define energies of conformations?
            ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
            return confEcalcMinimized;
        }
    }
}
