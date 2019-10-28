package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.pruning.AStarPruner;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseRigidGScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.StaticBiggestLowerboundDifferenceOrder;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.sharkstar.tools.SHARKStarEnsembleAnalyzer;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;
import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.printTree;

public class MultiSequenceSHARKStarBound implements PartitionFunction {

    private final EnergyMatrixCorrector energyMatrixCorrector = new EnergyMatrixCorrector(this);
    private Sequence precomputedSequence;
    protected double targetEpsilon = 1;
    public static final boolean debug = true;
    public boolean profileOutput = false;
    private Status status = null;

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
    public ArrayList<Integer> minList;
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

    private SHARKStarEnsembleAnalyzer ensembleAnalyzer;
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
    public final SimpleConfSpace confSpace;

    private BigDecimal precomputedUpperBound;
    private BigDecimal precomputedLowerBound;

    private List<MultiSequenceSHARKStarNode> precomputedFringe = new ArrayList<>();

    public static final int[] debugConf = new int[]{-1, 5, 3, 3, 8, 3};//4, -1, 8, 9};
    private static final int[] debugConf2 = new int[]{};//{8, 3, 4, 5, 8};
    private String cachePattern = "NOT_INITIALIZED";

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
        this.confSpace = confSpace;
        gscorerFactory = (emats) -> new PairwiseGScorer(emats);
        rigidgscorerFactory = (emats) -> new PairwiseRigidGScorer(emats);


        hscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, false);
        nhscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, true);
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
        double confUpperBound = partialConfUpperBound + nhscorerFactory.make(rigidEmat).calc(confIndex, rcs);
        rootNode.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperBound, precomputedSequence);
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.fullRCs = rcs;
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat));
        /* force init order */
        order.getNextPos(confIndex,fullRCs);
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
        ensembleAnalyzer = new SHARKStarEnsembleAnalyzer(minimizingEcalc, minimizingEmat);

        setParallelism(parallelism);

        // Recording pfunc starting bounds
        this.startLowerBound = rootNode.getLowerBound(precomputedSequence);
        this.startUpperBound = rootNode.getUpperBound(precomputedSequence);
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(), 0));
    }

    public MultiSequenceSHARKStarBound(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                       ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism,
                                       MultiSequenceSHARKStarBound precomputedFlex) {

        this(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, rcs, parallelism);
        processPrecomputedFlex(precomputedFlex);
    }

    private void processPrecomputedFlex(MultiSequenceSHARKStarBound precomputedFlex) {
        precomputedPfunc = precomputedFlex;
        precomputedRootNode = precomputedFlex.rootNode;
        this.precomputedSequence = precomputedFlex.confSpace.makeWildTypeSequence();
        precomputedUpperBound = precomputedRootNode.getUpperBound(precomputedSequence);
        precomputedLowerBound = precomputedRootNode.getLowerBound(precomputedSequence);
        updatePrecomputedConfTree();
        mergeCorrections(precomputedFlex.correctionMatrix, genConfSpaceMapping());

        // Fix order issues
        ConfIndex rootIndex = new ConfIndex(fullRCs.getNumPos());
        this.rootNode.getConfSearchNode().index(rootIndex);
        this.order.updateForPrecomputedOrder(precomputedFlex.order, rootIndex, this.fullRCs, genConfSpaceMapping());
    }

    /**
     * Takes partial minimizations from the precomputed correctionMatrix, maps them to the new confspace, and
     * stores them in this correctionMatrix
     */
    public void mergeCorrections(UpdatingEnergyMatrix precomputedCorrections, int[] confSpacePermutation){
        List<TupE> corrections = precomputedCorrections.getAllCorrections().stream()
                .map((tup) -> {
                    return tup.permute(confSpacePermutation);
                })
                .collect(Collectors.toList());
        if (corrections.size()!=0) {
            int TestNumCorrections = corrections.size();
            this.correctionMatrix.insertAll(corrections);
        }else
            System.out.println("No corrections to insert");
    }

    /**
     * Returns a wrapped pointer to this class, so that BBK* and MSK* can pretend they have single-sequence
     * partition functions.
     */
    public PartitionFunction getPartitionFunctionForSequence(Sequence seq) {
        SingleSequenceSHARKStarBound newBound = new SingleSequenceSHARKStarBound(this, seq, this);
        newBound.init(null, null, targetEpsilon);
        System.out.println("Creating new pfunc for sequence "+seq);
        System.out.println("Full RCs: "+fullRCs);
        System.out.println("Sequence RCs: "+newBound.seqRCs);
        computeFringeForSequence(newBound, this.rootNode);
        newBound.updateBound();
        if(newBound.getSequenceEpsilon() == 0)
            System.err.println("Perfectly bounded sequence? how?");
        rootNode.updateSubtreeBounds(seq);
        //printTree(seq, this.rootNode);
        //printTree("", null, confSpace, seq, this.rootNode);
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

    private void computeFringeForSequence(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarNode curNode) {
        Node confNode = curNode.getConfSearchNode();
        RCs rcs = bound.seqRCs;
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            confNode.index(context.index);

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

        System.out.println("Making flexible confspace bound...");

        ConfEnergyCalculator flexMinimizingConfECalc = FlexEmatMaker.makeMinimizeConfEcalc(flexConfSpace,
                this.parallelism);
        ConfEnergyCalculator rigidConfECalc = FlexEmatMaker.makeRigidConfEcalc(flexMinimizingConfECalc);
        EnergyMatrix flexMinimizingEmat = FlexEmatMaker.makeEmat(flexMinimizingConfECalc, "minimized", cachePattern+".flex");
        EnergyMatrix flexRigidEmat = FlexEmatMaker.makeEmat(rigidConfECalc, "rigid", cachePattern+".flex");
        UpdatingEnergyMatrix flexCorrection = new UpdatingEnergyMatrix(flexConfSpace, flexMinimizingEmat,
                flexMinimizingConfECalc);


        MultiSequenceSHARKStarBound precompFlex = new MultiSequenceSHARKStarBound(
                flexConfSpace, flexRigidEmat, flexMinimizingEmat,
                flexMinimizingConfECalc, flexRCs, this.parallelism);
        precompFlex.setCachePattern(cachePattern);
        precompFlex.initFlex(this.targetEpsilon, this.stabilityThreshold, flexCorrection);
        PartitionFunction flexBound =
                precompFlex.getPartitionFunctionForSequence(unassignedFlex);
        SingleSequenceSHARKStarBound bound = (SingleSequenceSHARKStarBound) flexBound;
        flexBound.compute();
        precompFlex.printEnsembleAnalysis();
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
        System.out.println("Tightening bound for "+sequenceBound.sequence);
        debugPrint("Num conformations: " + sequenceBound.numConformations);
        sequenceBound.updateBound();
        double lastEps = sequenceBound.getSequenceEpsilon();
        if(lastEps == 0)
            System.err.println("???!");

        int previousConfCount = workDone();

        if (!sequenceBound.nonZeroLower()) {
            runUntilNonZero(sequenceBound);
            sequenceBound.updateBound();
        }

        while (sequenceBound.getSequenceEpsilon() > targetEpsilon &&
                workDone() - previousConfCount < maxNumConfs
                && isStable(stabilityThreshold, sequenceBound.sequence)) {
            debugPrint("Tightening from epsilon of " + sequenceBound.getSequenceEpsilon());
            if (debug) {
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                debugHeap(sequenceBound.fringeNodes);
                debugHeap(sequenceBound.leafQueue);
                debugHeap(sequenceBound.internalQueue);
                //printTree(sequenceBound.sequence,rootNode);
            }
            tightenBoundInPhases(sequenceBound);
            debugPrint("Errorbound is now " + sequenceBound.getSequenceEpsilon());
            debugPrint("Bound reduction: "+(lastEps - sequenceBound.getSequenceEpsilon()));
            if (lastEps < sequenceBound.getSequenceEpsilon() && sequenceBound.getSequenceEpsilon() - lastEps > 0.01
                || sequenceBound.errors()) {
                System.err.println("Error. Bounds got looser.");
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                //printTree(sequenceBound.sequence,rootNode);
                System.exit(-1);
            }
            lastEps = sequenceBound.getSequenceEpsilon();
        }
        if (!isStable(stabilityThreshold, sequenceBound.sequence))
            sequenceBound.setStatus(Status.Unstable);
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
        throw new UnsupportedOperationException("Multisequence results are ill-defined.");
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
        if (queue.isEmpty()) {
            if(!bound.leafQueue.isEmpty()&&
                    MathTools.isGreaterThan(bound.leafQueue.peek().getLowerBound(bound.sequence),
                            BigDecimal.ONE))
                queue.add(bound.leafQueue.poll());
            else if(!bound.internalQueue.isEmpty() )
                queue.add(bound.internalQueue.poll());
            else
                queue.add(rootNode);
        }
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
                bound.getSequenceEpsilon(), bound.getValues().calcLowerBound(),
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
        //debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]", bound.getValues().calcLowerBound(),
        //        bound.getValues().calcUpperBound()));
        internalZ = ZSums[0];
        leafZ = ZSums[1];
        if(MathTools.isRelativelySame(internalZ, leafZ, PartitionFunction.decimalPrecision, 1e-3)
                && MathTools.isRelativelySame(leafZ, BigDecimal.ZERO, PartitionFunction.decimalPrecision, 1e-3)) {
            rootNode.updateSubtreeBounds(bound.sequence);
            printTree(bound.sequence, rootNode);
            System.out.println("This is a bad time.");
        }
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        if(!bound.internalQueue.isEmpty() &&
                MathTools.isLessThan(internalZ, bound.internalQueue.peek().getUpperBound(bound.sequence)))
            System.out.println("Should have used a node from the internal queue. How??");
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
        if (leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int) Math.floor(0.1 * leafTimeAverage / internalTimeAverage));
        while (!queue.isEmpty() && (bound.internalQueue.size() < maxNodes
                || (!bound.leafQueue.isEmpty() && MathTools.isGreaterThan(queue.peek().getErrorBound(bound.sequence),
                                                                            bound.leafQueue.peek().getErrorBound()))
                || bound.leafQueue.size() < maxMinimizations)) {
            MultiSequenceSHARKStarNode curNode = queue.poll();
            if(confMatch(debugConf, curNode.getConfSearchNode().assignments))
                System.out.println("Gotcha-populate");
            Node node = curNode.getConfSearchNode();
            ConfIndex index = new ConfIndex(fullRCs.getNumPos());

            node.index(index);
            double correctgscore = correctionMatrix.confE(node.assignments);
            double hscore = curNode.getConfLowerBound(bound.sequence) - node.getPartialConfLowerBound();
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
            if(confMatch(debugConf, curNode.getConfSearchNode().assignments))
                System.out.println("Gotcha-fillList");
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
        energyMatrixCorrector.processPreminimization(bound, minimizingEcalc);
        profilePrint("Preminimization time : " + loopWatch.getTime(2));
        double curEpsilon = bound.getSequenceEpsilon();
        //printTree(bound.sequence, rootNode);
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
            curNode.setBoundsFromConfLowerAndUpper(curNode.getConfLowerBound(seq), curNode.getConfUpperBound(seq), seq);
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
                    progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
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
                    progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
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
        System.out.println("Bounding "+startNode.toSeqString(bound.sequence));
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
                // Sometimes there are no good leaf nodes. Weird.
                if(nextNode != null) {
                    newNodes.remove(nextNode);
                    drillQueue.add(nextNode);
                }
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

    /**
     * Returns a correction matrix with full minimizations included
     */
    public UpdatingEnergyMatrix genCorrectionMatrix() {
        addFullMinimizationsToCorrectionMatrix();
        return this.correctionMatrix;
    }

    /**
     * Takes the full minimizations from this tree, insert them into the correction matrix
     */
    private void addFullMinimizationsToCorrectionMatrix(){
        captureSubtreeFullMinimizations(this.rootNode);
    }

    /**
     * Takes the full minimizations from this subtree, insert them into the correction matrix
     */
    private void captureSubtreeFullMinimizations(MultiSequenceSHARKStarNode subTreeRoot){
        if (subTreeRoot.getChildren(precomputedSequence) == null || subTreeRoot.getChildren(precomputedSequence).size() ==0){
            if (subTreeRoot.isMinimized(precomputedSequence)){
                RCTuple tuple = subTreeRoot.toTuple();
                double confEnergy = subTreeRoot.getConfLowerBound(precomputedSequence);
                double lowerbound = this.minimizingEmat.getInternalEnergy(tuple);
                if (lowerbound == confEnergy)
                    throw new ValueException("Minimized energies shouldn't equal lower bounds");
                double correction = confEnergy - lowerbound;
                this.correctionMatrix.setHigherOrder(tuple, correction);
            }
        }else{
            for (MultiSequenceSHARKStarNode node : subTreeRoot.getChildren(precomputedSequence)){
                captureSubtreeFullMinimizations(node);
            }
        }

    }

    public List<TupE> getCorrections() {
        return correctionMatrix.getAllCorrections().stream()
                .map((tup) -> {
                    return tup.permute(genConfSpaceMapping());
                })
                .collect(Collectors.toList());
    }

    public int getNumPartialMinimizations() {
        return numPartialMinimizations;
    }

    public MARKStarProgress getProgress() {
        return progress;
    }

    public EnergyMatrix getRigidEmat() {
        return rigidEmat;
    }

    public boolean isComputedCorrections() {
        return computedCorrections;
    }

    public EnergyMatrix getMinimizingEmat() {
        return minimizingEmat;
    }

    public ConfEnergyCalculator getMinimizingEcalc() {
        return minimizingEcalc;
    }

    public Set<String> getCorrectedTuples() {
        return correctedTuples;
    }

    public UpdatingEnergyMatrix getCorrectionMatrix() {
        return correctionMatrix;
    }

    public void setNumPartialMinimizations(int numPartialMinimizations) {
        this.numPartialMinimizations = numPartialMinimizations;
    }

    public void setComputedCorrections(boolean computedCorrections) {
        this.computedCorrections = computedCorrections;
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

                    // score the child node differentially against the parent node
                    if (child.getLevel() < RCs.getNumPos()) {
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double diff = Math.max(confCorrection,
                                context.partialConfLowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc));
                        double rigiddiff = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double hdiff = context.lowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double maxhdiff = context.upperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double checkNum = 0;
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
                        progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
                        result.lowerBound = confLowerBound;
                        result.upperBound = confUpperbound;
                    }
                    if (child.getLevel() == RCs.getNumPos()) {
                        double confRigid = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        //confRigid = confRigid - node.partialConfLowerbound + node.partialConfUpperBound;

                        child.computeNumConformations(RCs); // Shouldn't this always eval to 1, given that we are looking at leaf nodes?
                        double confLower = context.partialConfLowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                        double confCorrection = correctionMatrix.confE(child.assignments);
                        double lowerbound = Math.max(minimizingEmat.confE(child.assignments), confLower);

                        if (lowerbound < confCorrection) {
                            recordCorrection(lowerbound, confCorrection - lowerbound);
                        }
                        checkBounds(confCorrection, confRigid);
                        child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                        child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                        numConfsScored++;
                        progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
                        result.lowerBound = confCorrection;
                        result.upperBound = confRigid;
                    }
                    partialTime.stop();
                    loopPartialTime += partialTime.getTimeS();


                    return result;
                }

            }, (PartialConfNodeResult result) -> {
                assert (result.isValid());
                MultiSequenceSHARKStarNode newChild = curNode.makeChild(result.resultNode,
                        bound.sequence, result.lowerBound, result.upperBound);
                newChild.setBoundsFromConfLowerAndUpper(result.lowerBound,
                        result.upperBound, bound.sequence);
                //System.out.println("Created new child "+MultiSequenceSHARKStarNodeChild.toSeqString(bound.sequence));
                // collect the possible children
                if (newChild.getConfLowerBound(bound.sequence) < 0) {
                    children.add(newChild);

                }
                if (!newChild.isMinimized(bound.sequence)) {
                    newNodes.add(newChild);
                } else {
                    newChild.computeEpsilonErrorBounds(bound.sequence);
                    bound.addFinishedNode(newChild);
                }

                //curNode.updateSubtreeBounds(bound.sequence);
                //printTree(bound.sequence, curNode);
                curNode.markUpdated();
            });
        }
    }

    public static boolean confMatch(int[] assignments, int[] ints) {
        if(assignments.length != ints.length)
            return false;
        for(int i = 0; i < assignments.length; i++) {
           if(assignments[i] != ints[i])
               return false;
        }
        return true;
    }

    private boolean confSubset(int[] assignments, int[] ints) {
        if(assignments.length != ints.length)
            return false;
        for(int i = 0; i < assignments.length; i++) {
            if(assignments[i] != ints[i] && assignments[i] >= 0)
                return false;
        }
        return true;
    }


    protected void processFullConfNode(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> newNodes,
                                       MultiSequenceSHARKStarNode curNode, Node node) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        if(curNode.getConfLowerBound(bound.sequence) > 10 &&
                (!bound.fringeNodes.isEmpty() && bound.fringeNodes.peek().getConfLowerBound(bound.sequence) < 0
                || !bound.internalQueue.isEmpty() && bound.internalQueue.peek().getConfLowerBound(bound.sequence) < 0
                || !bound.leafQueue.isEmpty() && bound.leafQueue.peek().getConfLowerBound(bound.sequence) < 0)) {
            System.out.println("not processing high-energy conformation");
            newNodes.add(curNode);
            return;
        }
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
                        energyMatrixCorrector.computeEnergyCorrection(analysis, conf, context.ecalc, bound.getSequenceEpsilon());

                        double energy = analysis.epmol.energy;
                        double newConfUpper = energy;
                        double newConfLower = energy;
                        // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
                        double oldConfUpper = curNode.getConfUpperBound(bound.sequence);
                        double oldConfLower = curNode.getConfLowerBound(bound.sequence);
                        if (newConfUpper > oldConfUpper) {
                            System.err.println("Upper bounds got worse after minimization:" + newConfUpper
                                    + " > " + (oldConfUpper) + ". Rejecting minimized energy.");
                            System.err.println("Node info: " + curNode.toSeqString(bound.sequence));

                            newConfUpper = oldConfUpper;
                            newConfLower = oldConfUpper;
                        }
                        curNode.setBoundsFromConfLowerAndUpper(newConfLower, newConfUpper, bound.sequence);
                        double oldgscore = node.getPartialConfLowerBound();
                        node.setPartialConfLowerAndUpper(newConfLower, newConfUpper);
                        String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (curNode.getConfLowerBound(bound.sequence)) + "," + (curNode.getConfUpperBound(bound.sequence)) + "]";
                        debugPrint(out);
                        //ensembleAnalyzer.analyzeFullConf(analysis, conf);
                        curNode.markUpdated();
                        synchronized (this) {
                            if(precomputedSequence.equals(confSpace.makeUnassignedSequence()))
                                correctionMatrix.setHigherOrder(curNode.toTuple(),
                                        energy - minimizingEmat.confE(node.assignments));
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
                    progress.reportLeafNode(node.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
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
                    seqBound.getSequenceEpsilon(), stopwatch.getTime(2)));

        }
    }


    private void checkBounds(double lower, double upper) {
        if (upper < lower && upper - lower > 1e-5 && upper < 10)
            debugPrint("Bounds incorrect.");
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

    public interface ScorerFactory {
        AStarScorer make(EnergyMatrix emat);
    }

    public void setCachePattern(String pattern){
        this.cachePattern = pattern;
    }

    public void printEnsembleAnalysis() {
        ensembleAnalyzer.printStats();
    }


}
