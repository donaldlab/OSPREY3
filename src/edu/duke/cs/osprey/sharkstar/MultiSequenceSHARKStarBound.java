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
import edu.duke.cs.osprey.energy.BatchCorrectionMinimizer;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.StaticBiggestLowerboundDifferenceOrder;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode.Node;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics;
import edu.duke.cs.osprey.sharkstar.tools.SHARKStarEnsembleAnalyzer;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;
import jdk.nashorn.internal.runtime.regexp.joni.exception.ValueException;
import org.graalvm.compiler.core.common.type.ArithmeticOpTable;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.printTree;

public class MultiSequenceSHARKStarBound implements PartitionFunction {

    private final EnergyMatrixCorrector energyMatrixCorrector;
    private Sequence precomputedSequence;
    protected double targetEpsilon = 1;
    public static final boolean debug = false;
    public static final boolean batcher = true;
    public boolean profileOutput = false;
    private Status status = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;
    private int numConfsEnergiedThisLoop = 0;
    // max confs minimized, -1 means infinite.
    private int maxNumConfs = -1;

    // the number of full conformations scored OR energied
    private int numConfsScored = 0;

    protected int numInternalNodesProcessed = 0;

    private boolean printMinimizedConfs = false;
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
    public ObjectPool<ScoreContext> contexts;
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
    final BatchCorrectionMinimizer theBatcher;

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
    private MultiSequenceState state;


    private MultiSequenceSHARKStarBound precomputedPfunc;
    public MultiSequenceSHARKStarNode precomputedRootNode;
    public final SimpleConfSpace confSpace;

    private BigDecimal precomputedUpperBound;
    private BigDecimal precomputedLowerBound;

    private List<MultiSequenceSHARKStarNode> precomputedFringe = new ArrayList<>();

    public static final int[] debugConf = new int[]{};//6, 5, 15, -1, -1, 8, -1};
    private boolean internalQueueWasEmpty = false;
    private Map<Sequence, List<String>> scoreHistory = new HashMap<>();
    private String cachePattern = "NOT_INITIALIZED";

    public static final boolean writeTimes = true;
    private BufferedWriter writer;

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
        this.minimizingEmat = minimizingEmat;
        this.rigidEmat = rigidEmat;
        this.fullRCs = rcs;
        this.pruner = null;
        this.correctionMatrix = new UpdatingEnergyMatrix(confSpace, minimizingEmat);

        this.state = new MultiSequenceState();

        confIndex = new ConfIndex(rcs.getNumPos());
        Node rootConfNode = new Node(confSpace.positions.size(), 0, new MathTools.DoubleBounds());
        rootConfNode.index(confIndex);
        double partialConfLowerbound = gscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        double partialConfUpperBound = rigidgscorerFactory.make(rigidEmat).calc(confIndex, rcs);
        rootConfNode.computeNumConformations(rcs);

        rootConfNode.setBoundsFromConfLowerAndUpper(partialConfLowerbound, partialConfUpperBound);

        // Initialize residue ordering
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat));
        /* force init order */
        order.getNextPos(confIndex,fullRCs);

        // No precomputed sequence means the "precomputed" sequence is empty
        this.precomputedSequence = confSpace.makeUnassignedSequence();

        this.rootNode = MultiSequenceSHARKStarNode.makeRoot(rootConfNode, confSpace,
                confSpace.positions.get(order.getNextPos(confIndex, rcs)));
        double confLowerBound = partialConfLowerbound + hscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        double confUpperBound = partialConfUpperBound + nhscorerFactory.make(rigidEmat).calc(confIndex, rcs);

        rootNode.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, precomputedSequence, "(root initialization)");

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
        this.theBatcher = new BatchCorrectionMinimizer(minimizingConfEcalc, correctionMatrix, minimizingEmat);

        progress = new MARKStarProgress(fullRCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        ensembleAnalyzer = new SHARKStarEnsembleAnalyzer(minimizingEcalc, minimizingEmat);
        energyMatrixCorrector = new EnergyMatrixCorrector(this);

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

        /*
        // Compute the fringe in parallel
        List<MultiSequenceSHARKStarNode> scoredFringe = computeFringeForSequenceParallel(newBound.sequence, newBound.seqRCs);
        if(scoredFringe.size()==0)
            scoredFringe.add(this.rootNode);
        debugPrint(String.format("[Normal fringe # nodes, Parallel fringe # nodes] = [%d, %d]",newBound.fringeNodes.size(), scoredFringe.size()));
        newBound.fringeNodes.addAll(scoredFringe);
         */

        //rootNode.updateSubtreeBounds(seq);
        newBound.updateBound();
        if(newBound.getSequenceEpsilon() == 0)
            System.err.println("Perfectly bounded sequence? how?");
        //rootNode.updateSubtreeBounds(seq);
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
        node.makeNodeCompatibleWithConfSpace(permutation, size, this.fullRCs);
        node.setNewConfSpace(confSpace);
        if (node.hasChildren(precomputedSequence)) {
            for (MultiSequenceSHARKStarNode child : node.getChildren(precomputedSequence)) {
                updatePrecomputedNode(child, permutation, size);
            }
        }else {
            precomputedFringe.add(node);
        }
    }
    private List<MultiSequenceSHARKStarNode> computeFringeForSequenceParallel(Sequence seq, RCs seqRCs) {
        List<MultiSequenceSHARKStarNode> scoredSeqFringe = Collections.synchronizedList(new ArrayList<>());
        ConcurrentLinkedQueue<MultiSequenceSHARKStarNode> unscoredSeqFringe = new ConcurrentLinkedQueue<>(precomputedFringe);

        // Sometimes the precomputedFringe will be empty. It really shouldn't be, but for now just add the root Node if it is
        if(unscoredSeqFringe.isEmpty())
            unscoredSeqFringe.add(this.rootNode);

        while(!unscoredSeqFringe.isEmpty() || loopTasks.isWorking()){
            loopTasks.submit(
                ()->{
                    MultiSequenceSHARKStarNode node = unscoredSeqFringe.poll();
                    // If queue is empty, just return
                    if(node == null)
                        return null;
                    // Else, do things
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        ScoreContext context = checkout.get();

                        // Fix issue where nextDesignPosition can be null
                        MultiSequenceSHARKStarNode.Node confNode = node.getConfSearchNode();
                        confNode.index(context.index);
                        if (node.nextDesignPosition == null && node.level < confSpace.positions.size()) {
                            node.nextDesignPosition = confSpace.positions.get(order.getNextPos(context.index, seqRCs));
                        }

                        // Get the children
                        List<MultiSequenceSHARKStarNode> children = node.getOrMakeChildren(seq);
                        // If we are at a leaf, score the node
                        if(!children.isEmpty()) {
                            // if there are children, just add them to queue, since we only want the fringe
                            unscoredSeqFringe.addAll(children);
                        }else{
                            double confCorrection = correctionMatrix.confE(confNode.assignments);
                            double gscore = context.partialConfLowerBoundScorer.calc(context.index, seqRCs);
                            double hscore = context.lowerBoundScorer.calc(context.index, seqRCs);
                            double confLowerBound = confNode.getPartialConfLowerBound() + context.lowerBoundScorer.calc(context.index, seqRCs);
                            double confUpperBound = confNode.getPartialConfUpperBound() + context.upperBoundScorer.calc(context.index, seqRCs);
                            String historyString = String.format("%s: previous lower bound %f, g score %f, hscore %f, f score %f corrected score %f, from %s",
                                    confNode.confToString(), node.getConfLowerBound(seq), gscore, hscore, gscore + hscore, confCorrection, getStackTrace());
                            node.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, seq, historyString);

                            if (node.getChildren(null).isEmpty())
                                correctionMatrix.setHigherOrder(node.toTuple(), confNode.getPartialConfLowerBound()
                                        - minimizingEmat.confE(confNode.assignments));

                            //bound.fringeNodes.add(node);
                            scoredSeqFringe.add(node);
                        }
                    }
                    return null;
                },
                (ignored)->{}
            );
        }
        loopTasks.waitForFinish();
        return scoredSeqFringe;

    }

    private void computeFringeForSequence(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarNode curNode) {
        // Fix issue where next design position may not be defined
        RCs rcs = bound.seqRCs;
        Node confNode = curNode.getConfSearchNode();
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            confNode.index(context.index);
            if (curNode.nextDesignPosition == null && curNode.level < confSpace.positions.size()) {
                curNode.nextDesignPosition = confSpace.positions.get(order.getNextPos(context.index, rcs));
            }
        }
        // Get the children
        List<MultiSequenceSHARKStarNode> children = curNode.getOrMakeChildren(bound.sequence);
        // If we are at a leaf, score the node
        if(children.isEmpty()){
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                confNode.index(context.index);

                double confCorrection = correctionMatrix.confE(confNode.assignments);
                double gscore = context.partialConfLowerBoundScorer.calc(context.index, rcs);
                double hscore = context.lowerBoundScorer.calc(context.index, rcs);
                double confLowerBound = confNode.getPartialConfLowerBound() + context.lowerBoundScorer.calc(context.index, rcs);
                double confUpperBound = confNode.getPartialConfUpperBound() + context.upperBoundScorer.calc(context.index, rcs);
                String historyString = String.format("%s: previous lower bound %f, g score %f, hscore %f, f score %f corrected score %f, from %s",
                        confNode.confToString(), curNode.getConfLowerBound(bound.sequence), gscore, hscore, gscore + hscore, confCorrection, getStackTrace());
                curNode.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, bound.sequence, historyString);
                if (curNode.getChildren(null).isEmpty())
                    correctionMatrix.setHigherOrder(curNode.toTuple(), confNode.getPartialConfLowerBound()
                            - minimizingEmat.confE(confNode.assignments));
        /*
        if(bound.sequence.countMutations() < 1 && !curNode.hasChildren(bound.sequence) && curNode.getChildren(null).size() > 0) {
            System.out.println("Gotta be careful here.");
            printTree(bound.sequence, rootNode);
            computeFringeForSequence(bound, curNode);
        }
         */
            }

            bound.fringeNodes.add(curNode);
        }
        else
            for(MultiSequenceSHARKStarNode child: children)
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
        // Sometimes our designs don't have immutable residues on one side.
        if(flexConfSpace.positions.size() < 1)
            return null;
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

    /**
     * Try to apply corrections to node.
     */
    private boolean correctNodeOrFalse(MultiSequenceSHARKStarNode node, SingleSequenceSHARKStarBound bound) {
        boolean corrected = false;
        double confCorrection = correctionMatrix.confE(node.getConfSearchNode().assignments);
        double oldg = node.getConfSearchNode().getPartialConfLowerBound();
        double correctionDiff = confCorrection - oldg;

        if ( correctionDiff > 1e-5) {
            corrected = true;

            BigDecimal oldZUpperBound = node.getUpperBound(bound.sequence);
            double oldConfLowerBound = node.getConfLowerBound(bound.sequence);

            // update the node gscore
            node.getConfSearchNode().setPartialConfLowerAndUpper(confCorrection, node.getConfSearchNode().getPartialConfUpperBound());
            recordCorrection(oldg, correctionDiff);
            String historyString = String.format("%s: correction from %f to %f, from ",
                    node.getConfSearchNode().confToString(), oldg, confCorrection);

            // update the node total scores
            node.setBoundsFromConfLowerAndUpperWithHistory(oldConfLowerBound + correctionDiff, node.getConfUpperBound(bound.sequence), bound.sequence, historyString);
            node.markUpdated();
            System.out.println("Correcting " + node.toSeqString(bound.sequence) +" correction ="+(correctionDiff) );

            BigDecimal newUpperBound = node.getUpperBound(bound.sequence);
            BigDecimal diffUB = newUpperBound.subtract(oldZUpperBound, PartitionFunction.decimalPrecision);
            if(diffUB.compareTo(BigDecimal.ZERO) > 0){
                throw new RuntimeException();
            }
            if(node.getUpperBound(bound.sequence).compareTo(node.getLowerBound(bound.sequence)) < 0){
                System.err.println(String.format("Insane bounds: UB (%1.9e) < LB (%1.9e)",
                        node.getUpperBound(bound.sequence),
                        node.getLowerBound(bound.sequence)
                        ));
                throw new RuntimeException();
            }
            synchronized(bound.state){
                bound.state.upperBound = bound.state.upperBound.add(diffUB, PartitionFunction.decimalPrecision);
            }
        }
        return corrected;
    }

    @Override
    public void compute(int maxNumConfs) {
        throw new UnsupportedOperationException("Do not try to run Multisequence SHARK* bounds directly. Call " +
                "getPartitionFunctionForSequence() and use the generated single sequence bound.");
    }

    private static enum Step {
        None,
        Expand,
        ExpandInBatches,
        Energy,
        Partial,
    }

    private static class ExpansionResult{
        List <MultiSequenceSHARKStarNode> newNodes;
        BigDecimal deltaLB;
        BigDecimal deltaUB;
        double timeS = 0;
        int numExpanded = 0;
    }

    public void computeForSequenceParallel(int maxNumConfs, SingleSequenceSHARKStarBound sequenceBound){
        System.out.println("Tightening bound for "+sequenceBound.sequence);

        // Test to make sure that we have a valid pfunc
        sequenceBound.updateBound();
        double lastEps = sequenceBound.getSequenceEpsilon();
        if(lastEps == 0)
            System.err.println("???!");

        int previousConfCount = workDone();

        // Run to populate the queues a little bit
        if (!sequenceBound.nonZeroLower()) {
            runUntilNonZero(sequenceBound);
            sequenceBound.updateBound();
        }

        // Initialize state
        sequenceBound.state.upperBound = sequenceBound.getUpperBound();
        sequenceBound.state.lowerBound = sequenceBound.getLowerBound();
        double curEps = sequenceBound.state.calcDelta();

        Step step = Step.None;

        /*
        Continue looping if 1) there are nodes in the fringe queue or
        2) we are running tasks that will result in nodes being added to the queue
         */
        this.numConfsEnergiedThisLoop = 0;

        while( (!sequenceBound.fringeNodes.isEmpty() || loopTasks.isExpecting())){

            synchronized(sequenceBound) {
                //TODO: apply partial minimizations
                sequenceBound.updateBound();

                // Early termination
                    curEps = sequenceBound.state.calcDelta();

                System.out.println(String.format("Epsilon: %.9f, Bounds:[%1.3e, %1.3e]",
                        //sequenceBound.getSequenceEpsilon(),
                        //sequenceBound.getLowerBound(),
                        //sequenceBound.getUpperBound()
                        curEps,
                        sequenceBound.state.lowerBound,
                        sequenceBound.state.upperBound
                ));

                if (curEps < targetEpsilon ||
                        numConfsEnergiedThisLoop >= maxNumConfs ||
                        !isStable(stabilityThreshold, sequenceBound)
                ){
                    if(workDone() - previousConfCount >= maxNumConfs)
                        System.out.println("Exiting loop because of work");

                    if(!isStable(stabilityThreshold, sequenceBound)) {
                        System.out.println("Exiting loop due to stablity, thresh: " + stabilityThreshold);
                        sequenceBound.setStatus(Status.Unstable);
                    }
                    break;
                }
            }


            // Make lists
            List<MultiSequenceSHARKStarNode> internalNodes = new ArrayList<>();
            List<MultiSequenceSHARKStarNode> leftoverLeaves = new ArrayList<>();
            MultiSequenceSHARKStarNode loosestLeaf = null;
            MultiSequenceSHARKStarNode loosestInternal = null;

            // Make stopwatches
            Stopwatch loopWatch = new Stopwatch();
            loopWatch.start();
            Stopwatch internalTime = new Stopwatch();

            int numNodes = 0;

            // How many internal nodes can we score in the time it takes to do a minimization?
            int maxInternalNodes = 1;
            //TODO: Figure out why setting the above to 1000 makes the main thread take all the resources
            synchronized(sequenceBound){
                double leafTimeAvg = sequenceBound.state.totalTimeEnergy / sequenceBound.state.numEnergiedConfs;
                double internalTimeAvg = sequenceBound.state.totalTimeExpansion / sequenceBound.state.numExpansions;
                if (leafTimeAverage > 0)
                    maxInternalNodes = Math.min(maxInternalNodes, (int) Math.floor(0.1 * leafTimeAvg / internalTimeAvg));
                maxInternalNodes = Math.min(maxInternalNodes,
                        (int) (sequenceBound.fringeNodes.size() / 2*loopTasks.getParallelism()));
                maxInternalNodes = Math.max(maxInternalNodes,1);
            }

            boolean computePartials = false;
            synchronized(theBatcher) {
                //Figure out what step to make and get nodes
                computePartials = theBatcher.canBatch();

                if(computePartials)
                    theBatcher.makeBatch();
                    step = Step.Partial;
            }
            synchronized(sequenceBound.fringeNodes) {
                if(!computePartials && sequenceBound.fringeNodes.size() > 0) {
                    MultiSequenceSHARKStarNode node = null;

                    boolean foundUncorrectedNode = false;
                    //while (!foundUncorrectedNode && internalNodes.size() < maxInternalNodes) {
                    while (!sequenceBound.fringeNodes.isEmpty() && internalNodes.size() < maxInternalNodes) {
                        node = sequenceBound.fringeNodes.poll();
                        // Try to correct the node
                        //foundUncorrectedNode = !correctNodeOrFalse(node, sequenceBound);
                        //foundUncorrectedNode = true;
                        // If we corrected the node, put it back and try again
                        //if (!foundUncorrectedNode) {
                            //sequenceBound.fringeNodes.add(node);
                            // If we didn't correct the node
                        //} else {
                            // if it's a leaf
                            if (node.getConfSearchNode().getLevel() >= fullRCs.getNumPos()) {
                                // if we don't already have a leaf
                                if (loosestLeaf == null) {
                                    loosestLeaf = node;
                                    // if we do already have the leaf, put it back later
                                } else {
                                    leftoverLeaves.add(node);
                                }
                                // if it's not a leaf, add it to the internal list
                            } else {
                                internalNodes.add(node);
                            }

                        //}
                    }
                    sequenceBound.fringeNodes.addAll(leftoverLeaves);
                    // now, decide which step to take
                    BigDecimal leafSum = BigDecimal.ZERO;
                    for (int i = 0; i < internalNodes.size(); i++) {
                        leafSum = leafSum.add(internalNodes.get(i).getErrorBound(sequenceBound.sequence), PartitionFunction.decimalPrecision);
                    }

                    if (loosestLeaf == null) {
                        step = Step.ExpandInBatches;
                    }else if(leafSum.compareTo(loosestLeaf.getErrorBound(sequenceBound.sequence)) > 0){
                        step = Step.ExpandInBatches;
                        sequenceBound.fringeNodes.add(loosestLeaf);
                    }else{
                        step = Step.Energy;
                        sequenceBound.fringeNodes.addAll(internalNodes);
                    }
                }
            }


            // Take steps
            switch(step) {
                case Energy:{
                    MultiSequenceSHARKStarNode toMinimize = loosestLeaf;

                    loopTasks.submit( () -> {
                        return minimizeNode(toMinimize, sequenceBound.sequence);
                            },
                            (result) -> {
                                synchronized(sequenceBound.state) {
                                    // Update partition function values
                                    sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
                                    sequenceBound.state.lowerBound = sequenceBound.state.lowerBound.add(result.deltaLB, PartitionFunction.decimalPrecision);
                                    if(result.didMinimize){
                                        sequenceBound.state.numEnergiedConfs++;
                                        sequenceBound.state.totalTimeEnergy+=result.timeS;
                                        sequenceBound.state.numRoundsEnergy++;
                                    }
                                }

                                synchronized (this) { // don't race the main thread
                                    if(result.didMinimize) {
                                        if (precomputedSequence.equals(confSpace.makeUnassignedSequence()))
                                            correctionMatrix.setHigherOrder(result.minimizedNode.toTuple(),
                                                    result.energy - minimizingEmat.confE(result.minimizedNode.getConfSearchNode().assignments));
                                        numConfsEnergied++;
                                        this.numConfsEnergiedThisLoop++;
                                        minList.set(result.conf.getAssignments().length - 1, minList.get(result.conf.getAssignments().length - 1) + 1);

                                        this.state.numEnergiedConfs++;
                                        this.state.totalTimeEnergy += result.timeS;
                                        this.state.numRoundsEnergy++;

                                        //recordReduction(oldConfLower, oldConfUpper, energy);
                                        //printMinimizationOutput(node, newConfLower, oldgscore, bound);
                                        sequenceBound.addFinishedNode(result.minimizedNode);

                                        progress.reportLeafNode(result.minimizedNode.getConfSearchNode().getPartialConfLowerBound(), sequenceBound.fringeNodes.size(), sequenceBound.getSequenceEpsilon());
                                        leafTimeAverage = result.timeS;
                                        System.out.println("Processed 1 leaf in " + result.timeS + " seconds.");
                                    }else{
                                        synchronized(sequenceBound.fringeNodes){
                                            sequenceBound.fringeNodes.add(result.minimizedNode);

                                        }
                                    }
                                }
                            }
                    );
                    break;
                }
                case Expand: {
                    loopTasks.submitExpecting(
                            () -> {
                                ExpansionResult result = new ExpansionResult();
                                //result.newNodes = Collections.synchronizedList(new ArrayList<>());

                                internalTime.reset();
                                internalTime.start();

                                BigDecimal startLB = BigDecimal.ZERO;
                                BigDecimal startUB = BigDecimal.ZERO;
                                BigDecimal endLB = BigDecimal.ZERO;
                                BigDecimal endUB = BigDecimal.ZERO;

                                for(int j=0; j< 100; j++) {

                                    MultiSequenceSHARKStarNode node;
                                    synchronized (sequenceBound.fringeNodes) {
                                        if (sequenceBound.fringeNodes.isEmpty()) {
                                            break;
                                        }
                                        node = sequenceBound.fringeNodes.poll();

                                        // Make sure we have an internal node
                                        List<MultiSequenceSHARKStarNode> tempList = new ArrayList<>();
                                        while (node.getConfSearchNode().getLevel() == fullRCs.getNumPos()) {
                                            tempList.add(node);
                                            node = sequenceBound.fringeNodes.poll();
                                        }
                                        sequenceBound.fringeNodes.addAll(tempList);
                                    }

                                    List<MultiSequenceSHARKStarNode> newNodes = new ArrayList<>();

                                    startLB = startLB.add(node.getLowerBound(sequenceBound.sequence), PartitionFunction.decimalPrecision);
                                    startUB = startUB.add(node.getUpperBound(sequenceBound.sequence), PartitionFunction.decimalPrecision);


                                    node.checkDescendents(sequenceBound.sequence);
                                    if (!MathTools.isGreaterThan(node.getLowerBound(sequenceBound.sequence), BigDecimal.ONE) &&
                                            MathTools.isGreaterThan(
                                                    MathTools.bigDivide(node.getUpperBound(sequenceBound.sequence), sequenceBound.state.upperBound,
                                                            PartitionFunction.decimalPrecision),
                                                    new BigDecimal(1 - targetEpsilon))
                                    ) {
                                        boundLowestBoundConfUnderNode(sequenceBound, node, newNodes);
                                    } else {
                                        processPartialConfNode(sequenceBound, newNodes, node, node.getConfSearchNode());
                                    }
                                    node.markUpdated();

                                    endLB = endLB.add(newNodes.stream()
                                            .map(n -> n.getLowerBound(sequenceBound.sequence))
                                            .reduce(BigDecimal.ZERO, (a, b) -> a.add(b, PartitionFunction.decimalPrecision)), PartitionFunction.decimalPrecision);
                                    endUB = endUB.add(newNodes.stream()
                                            .map(n -> n.getUpperBound(sequenceBound.sequence))
                                            .reduce(BigDecimal.ZERO, (a, b) -> a.add(b, PartitionFunction.decimalPrecision)), PartitionFunction.decimalPrecision);

                                    synchronized(sequenceBound.fringeNodes){
                                        sequenceBound.fringeNodes.addAll(newNodes);
                                    }
                                }

                            result.deltaLB = endLB.subtract(startLB, PartitionFunction.decimalPrecision);
                            result.deltaUB = endUB.subtract(startUB, PartitionFunction.decimalPrecision);

                            internalTime.stop();

                            return result;
                        },
                        (result) -> {
                            synchronized(sequenceBound.state) {
                                // Update partition function values
                                sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
                                sequenceBound.state.lowerBound = sequenceBound.state.lowerBound.add(result.deltaLB, PartitionFunction.decimalPrecision);
                            }

                            synchronized (this) {
                                System.out.println("Got some internal nodes.");
                                internalTimeSum = internalTime.getTimeS();
                                internalTimeAverage = internalTimeSum / Math.max(1, 1);
                                debugPrint("Internal node time :" + internalTimeSum + ", average " + internalTimeAverage);
                                numInternalNodesProcessed += internalNodes.size();
                                synchronized (sequenceBound) {
                                    //sequenceBound.fringeNodes.addAll(result.newNodes);

                                }
                            }
                        }
                            );
                    break;
                }
                case ExpandInBatches: {
                    List<MultiSequenceSHARKStarNode> toExpand = internalNodes;
                    numNodes = toExpand.size();
                    System.out.println("Processing " + numNodes + " internal nodes...");

                    loopTasks.submitExpecting(
                            () -> {
                                ExpansionResult result = new ExpansionResult();
                                internalTime.reset();
                                internalTime.start();

                                BigDecimal startLB = BigDecimal.ZERO;
                                BigDecimal startUB = BigDecimal.ZERO;

                                result.newNodes = Collections.synchronizedList(new ArrayList<>());
                                for (MultiSequenceSHARKStarNode internalNode : toExpand) {

                                    startLB = startLB.add(internalNode.getLowerBound(sequenceBound.sequence), PartitionFunction.decimalPrecision);
                                    startUB = startUB.add(internalNode.getUpperBound(sequenceBound.sequence), PartitionFunction.decimalPrecision);

                                    internalNode.checkDescendents(sequenceBound.sequence);
                                    if (!MathTools.isGreaterThan(internalNode.getLowerBound(sequenceBound.sequence), BigDecimal.ONE) &&
                                            (sequenceBound.state.upperBound.compareTo(BigDecimal.ONE) < 0 ||
                                            MathTools.isGreaterThan(
                                                    MathTools.bigDivide(internalNode.getUpperBound(sequenceBound.sequence), sequenceBound.state.upperBound,
                                                            PartitionFunction.decimalPrecision),
                                                    new BigDecimal(1 - targetEpsilon)))
                                    ) {
                                        boundLowestBoundConfUnderNode(sequenceBound, internalNode, result.newNodes);
                                    } else {
                                        processPartialConfNode(sequenceBound, result.newNodes, internalNode, internalNode.getConfSearchNode());
                                    }
                                    internalNode.markUpdated();
                                }

                                BigDecimal endLB = result.newNodes.stream()
                                        .map( n -> n.getLowerBound(sequenceBound.sequence))
                                        .reduce(BigDecimal.ZERO, (a,b) -> a.add(b, PartitionFunction.decimalPrecision));
                                BigDecimal endUB = result.newNodes.stream()
                                        .map( n -> n.getUpperBound(sequenceBound.sequence))
                                        .reduce(BigDecimal.ZERO, (a,b) -> a.add(b, PartitionFunction.decimalPrecision));

                                result.deltaLB = endLB.subtract(startLB, PartitionFunction.decimalPrecision);
                                result.deltaUB = endUB.subtract(startUB, PartitionFunction.decimalPrecision);

                                internalTime.stop();
                                result.timeS = internalTime.getTimeS();
                                result.numExpanded = toExpand.size();
                                return result;
                            },
                            (result) -> {
                                synchronized(sequenceBound.state) {
                                    // Update partition function values
                                    sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
                                    sequenceBound.state.lowerBound = sequenceBound.state.lowerBound.add(result.deltaLB, PartitionFunction.decimalPrecision);
                                    sequenceBound.state.numExpansions += result.numExpanded;
                                    sequenceBound.state.totalTimeExpansion += result.timeS;
                                    sequenceBound.state.numRoundsExpand++;
                                }

                                synchronized (this) {
                                    this.state.numExpansions += result.numExpanded;
                                    this.state.totalTimeExpansion += result.timeS;
                                    this.state.numRoundsExpand++;
                                    System.out.println("Got " + result.newNodes.size() +" internal nodes.");
                                    internalTimeSum = internalTime.getTimeS();
                                    internalTimeAverage = internalTimeSum / Math.max(1, toExpand.size());
                                    debugPrint("Internal node time :" + internalTimeSum + ", average " + internalTimeAverage);
                                    numInternalNodesProcessed += internalNodes.size();
                                    synchronized (sequenceBound.fringeNodes) {
                                        sequenceBound.fringeNodes.addAll(result.newNodes);

                                    }
                                }
                            }
                    );

                    break;
                }
                case Partial: {
                    System.out.println("Computing partial mins");
                    BatchCorrectionMinimizer.Batch batch = theBatcher.acquireBatch();
                    loopTasks.submit(
                            () -> {
                                PartialMinimizationResult result = new PartialMinimizationResult();
                                Stopwatch partialMinTime = new Stopwatch().start();

                                // calculate all the fragment energies
                                Map<RCTuple, EnergyCalculator.EnergiedParametricMolecule> confs = new HashMap<>();
                                for (RCTuple frag : batch.fragments) {

                                    double energy;

                                    // are there any RCs are from two different backbone states that can't connect?
                                    if (theBatcher.isParametricallyIncompatible(frag)) {

                                        // yup, give this frag an infinite energy so we never choose it
                                        energy = Double.POSITIVE_INFINITY;

                                    } else {

                                        // nope, calculate the usual fragment energy
                                        confs.put(frag, theBatcher.confEcalc.calcEnergy(frag));
                                    }
                                }
                                partialMinTime.stop();

                                result.confs = confs;
                                result.numPartials = confs.size();
                                result.timeS = partialMinTime.getTimeS();
                                return result;
                            },
                            (PartialMinimizationResult result) -> {
                                // update the energy matrix
                                for(RCTuple tuple : result.confs.keySet()) {
                                    double lowerbound = theBatcher.minimizingEnergyMatrix.getInternalEnergy(tuple);
                                    double tupleEnergy = result.confs.get(tuple).energy;
                                    if (tupleEnergy - lowerbound > 0) {
                                        double correction = tupleEnergy - lowerbound;
                                        correctionMatrix.setHigherOrder(tuple, correction);
                                    } else
                                        System.err.println("Negative correction for " + tuple.stringListing());
                                }

                                // Record stats
                                synchronized(sequenceBound.state){
                                    sequenceBound.state.numPartialMinimizations += result.numPartials;
                                    sequenceBound.state.totalTimePartialMin += result.timeS;
                                    sequenceBound.state.numRoundsPartialMin++;
                                }
                                synchronized(this){
                                    this.state.numPartialMinimizations += result.numPartials;
                                    this.state.totalTimePartialMin += result.timeS;
                                    this.state.numRoundsPartialMin++;
                                }
                            }
                    );
                    break;
                }
                case None: {
                }

            }

            synchronized (this){
                /*
                sequenceBound.updateBound();
                loopWatch.stop();
                double loopTime = loopWatch.getTimeS();
                //profilePrint("Processed " + numNodes + " this loop, spawning " + newNodes.size() + " in " + loopTime + ", " + stopwatch.getTime() + " so far");
                loopWatch.reset();
                loopWatch.start();
                //energyMatrixCorrector.processPreminimization(bound, minimizingEcalc);
                //profilePrint("Preminimization time : " + loopWatch.getTime(2));
                double curEpsilon = sequenceBound.getSequenceEpsilon();
                //printTree(bound.sequence, rootNode);
                loopWatch.stop();
                //cleanupTime = loopWatch.getTimeS();
                //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);
                System.out.println(String.format("Loop complete. Bounds are now [%12.6e,%12.6e]", sequenceBound.getValues().calcLowerBound(),
                        sequenceBound.getValues().calcUpperBound()));
                //loopCleanup(sequenceBound, newNodes, loopWatch, numNodes);

                 */
            }

        }
        System.out.println("exited loop");
        System.out.println(sequenceBound.fringeNodes.size());
        loopTasks.waitForFinish();
        System.out.println("Finished tasks");
        System.out.println(sequenceBound.fringeNodes.size());
        sequenceBound.updateBound();
        curEps = sequenceBound.state.calcDelta();
        System.out.println(String.format("Tracking Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                curEps,
                sequenceBound.state.lowerBound,
                sequenceBound.state.upperBound
        ));
        System.out.println(String.format("Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                sequenceBound.getSequenceEpsilon(),
                sequenceBound.getLowerBound(),
                sequenceBound.getUpperBound()
        ));
        System.out.println(String.format("Minimized %d nodes, %d nodes in fringe.",
                sequenceBound.finishedNodes.size(),
                sequenceBound.fringeNodes.size()));

        System.out.println(String.format("--- Minimizations --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
                this.state.numEnergiedConfs,
                this.state.totalTimeEnergy / this.state.numEnergiedConfs,
                this.state.totalTimeEnergy / this.state.numRoundsEnergy
                ));
        System.out.println(String.format("--- Expansions --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
                this.state.numExpansions,
                this.state.totalTimeExpansion / this.state.numExpansions,
                this.state.totalTimeExpansion / this.state.numRoundsExpand
        ));
        System.out.println(String.format("--- Partials --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
                this.state.numPartialMinimizations,
                this.state.totalTimePartialMin / this.state.numPartialMinimizations,
                this.state.totalTimePartialMin / this.state.numRoundsPartialMin
        ));

    }

    private static class MinimizationResult{
        boolean didMinimize; // did we minimize, or just correct?
        MultiSequenceSHARKStarNode minimizedNode; // To put back in the fringe queue
        double energy;
        ConfAnalyzer.ConfAnalysis analysis; // For the batch minimizer
        ConfSearch.ScoredConf conf; // For the batch minimizer
        double timeS; // Minimization time in seconds
        BigDecimal deltaLB; // Lower bound change
        BigDecimal deltaUB; // Upper bound change
    }

    private static class PartialMinimizationResult{
        Map<RCTuple, EnergyCalculator.EnergiedParametricMolecule> confs;
        int numPartials;
        double timeS;
    }

    private MinimizationResult minimizeNode(MultiSequenceSHARKStarNode node, Sequence sequence) {
        MinimizationResult result = new MinimizationResult();
        Stopwatch minTimer = new Stopwatch().start();

        BigDecimal startingLB = node.getLowerBound(sequence);
        BigDecimal startingUB = node.getUpperBound(sequence);

        // First, try to correct the node instead of minimize it
        boolean corrected = false;
        double confCorrection = correctionMatrix.confE(node.getConfSearchNode().assignments);
        double oldg = node.getConfSearchNode().getPartialConfLowerBound();
        double correctionDiff = confCorrection - oldg;

        if ( correctionDiff > 1e-5) {
            corrected = true;

            BigDecimal oldZUpperBound = node.getUpperBound(sequence);
            double oldConfLowerBound = node.getConfLowerBound(sequence);

            // update the node gscore
            node.getConfSearchNode().setPartialConfLowerAndUpper(confCorrection, node.getConfSearchNode().getPartialConfUpperBound());
            String historyString = String.format("%s: correction from %f to %f, from ",
                    node.getConfSearchNode().confToString(), oldg, confCorrection);

            // update the node total scores
            node.setBoundsFromConfLowerAndUpperWithHistory(oldConfLowerBound + correctionDiff, node.getConfUpperBound(sequence), sequence, historyString);
            node.markUpdated();
            //System.out.println("Correcting " + node.toSeqString(sequence) + " correction =" + (correctionDiff));

            BigDecimal newUpperBound = node.getUpperBound(sequence);
            BigDecimal diffUB = newUpperBound.subtract(oldZUpperBound, PartitionFunction.decimalPrecision);
            if (diffUB.compareTo(BigDecimal.ZERO) > 0) {
                throw new RuntimeException();
            }
            if (node.getUpperBound(sequence).compareTo(node.getLowerBound(sequence)) < 0) {
                System.err.println(String.format("Insane bounds: UB (%1.9e) < LB (%1.9e)",
                        node.getUpperBound(sequence),
                        node.getLowerBound(sequence)
                ));
                throw new RuntimeException();
            }
            minTimer.stop();

            result.minimizedNode = node;
            result.didMinimize = false;
            result.timeS = minTimer.getTimeS();
            result.deltaLB = BigDecimal.ZERO;
            result.deltaUB = diffUB;

        }else{

            // now, try to minimize it if the correction failed
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                node.index(context.index);

                //System.out.println("Minmized "+curNode.toSeqString(bound.sequence));
                ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.getConfSearchNode().assignments, node.getConfLowerBound(sequence));
                ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);

                if (batcher) {
                    energyMatrixCorrector.scheduleEnergyCorrection(analysis, conf,
                            theBatcher);
                } else {
                    //computeEnergyCorrectionWithoutBatcher(analysis, conf, context.ecalc);
                }


                double energy = analysis.epmol.energy;
                double newConfUpper = energy;
                double newConfLower = energy;
                // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
                double oldConfUpper = node.getConfUpperBound(sequence);
                double oldConfLower = node.getConfLowerBound(sequence);
                if (newConfUpper > oldConfUpper) {
                    System.err.println("Upper bounds got worse after minimization:" + newConfUpper
                            + " > " + (oldConfUpper) + ". Rejecting minimized energy.");
                    System.err.println("Node info: " + node.toSeqString(sequence));

                    newConfUpper = oldConfUpper;
                    newConfLower = oldConfUpper;
                }

                String historyString = "";
                if (debug) {
                    historyString = String.format("minimimized %s to %s from %s",
                            node.getConfSearchNode().confToString(), energy, getStackTrace());
                }

                node.setBoundsFromConfLowerAndUpperWithHistory(newConfLower, newConfUpper, sequence, historyString);
                double oldgscore = node.getConfSearchNode().getPartialConfLowerBound();
                node.getConfSearchNode().setPartialConfLowerAndUpper(newConfLower, newConfUpper);
                String out = "Energy = " + String.format("%6.3e", energy) + ", [" + (node.getConfLowerBound(sequence)) + "," + (node.getConfUpperBound(sequence)) + "]";
                debugPrint(out);
                node.markUpdated();

                BigDecimal endLB = node.getLowerBound(sequence);
                BigDecimal endUB = node.getUpperBound(sequence);

                minTimer.stop();

                result.didMinimize = true;
                result.minimizedNode = node;
                result.energy = energy;
                result.analysis = analysis;
                result.conf = conf;
                result.timeS = minTimer.getTimeS();
                result.deltaLB = endLB.subtract(startingLB, PartitionFunction.decimalPrecision);
                result.deltaUB = endUB.subtract(startingUB, PartitionFunction.decimalPrecision);

            }

        }
        return result;
    }


    public void computeForSequence(int maxNumConfs, SingleSequenceSHARKStarBound sequenceBound) {
        if(writeTimes) {
            try {
                writer = new BufferedWriter(new FileWriter(stateName.concat("_shark.debug")));
                writer.write("popQueues time, internal time, internal nodes, leaf time, leaf nodes, cleanup time, total time, epsilon change\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

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
                && isStable(stabilityThreshold, sequenceBound)) {
            debugPrint("Tightening from epsilon of " + sequenceBound.getSequenceEpsilon());
            if (debug) {
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                debugHeap(sequenceBound.fringeNodes);
                debugHeap(sequenceBound.leafQueue);
                debugHeap(sequenceBound.internalQueue);
                internalQueueWasEmpty = sequenceBound.internalQueue.isEmpty();
                rootNode.debugTree(sequenceBound.sequence);
            }
            Stopwatch loopTimer = new Stopwatch();
            loopTimer.start();
            tightenBoundInPhases(sequenceBound);
            loopTimer.stop();
            debugPrint("Errorbound is now " + sequenceBound.getSequenceEpsilon());
            double delEps = lastEps - sequenceBound.getSequenceEpsilon();
            debugPrint("Bound reduction: "+delEps);
            if (lastEps < sequenceBound.getSequenceEpsilon() && sequenceBound.getSequenceEpsilon() - lastEps > 0.01
                || sequenceBound.errors()) {
                System.err.println("Error. Bounds got looser.");
                rootNode.updateSubtreeBounds(sequenceBound.sequence);
                //rootNode.debugTree(sequenceBound.sequence);
                System.exit(-1);
            }
            lastEps = sequenceBound.getSequenceEpsilon();

            if(writeTimes) {
                try {
                    writer.write(String.format(", %f, %.10f\n", loopTimer.getTimeS(), delEps));
                    writer.flush();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        loopTasks.waitForFinish();
        if (!isStable(stabilityThreshold, sequenceBound))
            sequenceBound.setStatus(Status.Unstable);
        BigDecimal averageReduction = BigDecimal.ZERO;
        int totalMinimizations = numConfsEnergied + numPartialMinimizations;
        if (totalMinimizations > 0)
            averageReduction = cumulativeZCorrection
                    .divide(new BigDecimal(totalMinimizations), new MathContext(BigDecimal.ROUND_HALF_UP));
        debugPrint(String.format("Average Z reduction per minimization: %12.6e", averageReduction));

        if(writeTimes) {
            try {
                writer.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
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
        // Use the threadpool that was created for emat calculation
        loopTasks = minimizingEcalc.tasks;
        /*
        if (loopTasks == null)
            loopTasks = parallelism.makeTaskExecutor(null);
        contexts.allocate(parallelism.getParallelism());

         */
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
            if(!bound.leafQueue.isEmpty())
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
        //rootNode.computeEpsilonErrorBounds(bound.sequence);
        nonZeroLower = true;
    }

    protected void tightenBoundInPhases(SingleSequenceSHARKStarBound bound) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        assert (!bound.isEmpty());
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
        Stopwatch popQueuesTimer = new Stopwatch();
        Stopwatch internalTime = new Stopwatch();
        Stopwatch leafTime = new Stopwatch();
        Stopwatch cleanupTime = new Stopwatch();
        int numInternals = 0;
        int numLeaves = 0;
        double leafTimeSum = 0;
        double internalTimeSum = 0;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};

        popQueuesTimer.start();
        populateQueues(bound, internalNodes, leafNodes, ZSums);
        //bound.updateBound();
        //debugPrint(String.format("After corrections, bounds are now [%12.6e,%12.6e]", bound.getValues().calcLowerBound(),
        //        bound.getValues().calcUpperBound()));
        internalZ = ZSums[0];
        leafZ = ZSums[1];
        if(leafNodes.isEmpty() && internalNodes.isEmpty()) {
            System.out.println("Nothing was populated?");
            populateQueues(bound, internalNodes, leafNodes, ZSums);
        }
        if(MathTools.isRelativelySame(internalZ, leafZ, PartitionFunction.decimalPrecision, 1e-3)
                && MathTools.isRelativelySame(leafZ, BigDecimal.ZERO, PartitionFunction.decimalPrecision, 1e-3)) {
            rootNode.updateSubtreeBounds(bound.sequence);
            //printTree(bound.sequence, rootNode);
            System.out.println("This is a bad time.");
            populateQueues(bound, new ArrayList<>(), new ArrayList<>(), ZSums);
        }
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        popQueuesTimer.stop();

        if(!bound.internalQueue.isEmpty() &&
                MathTools.isLessThan(internalZ, bound.internalQueue.peek().getUpperBound(bound.sequence)))
            System.out.println("Should have used a node from the internal queue. How??");
        if (MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            numLeaves = numNodes;
            System.out.println("Processing " + numNodes + " leaf nodes...");
            leafTime.reset();
            leafTime.start();
            /*
            if(bound.maxMinimizations < parallelism.numThreads)
                bound.maxMinimizations++;
             */
            for (MultiSequenceSHARKStarNode leafNode : leafNodes) {
                //debugPrint("Processing Node: " + leafNode.toSeqString(bound.sequence));
                processFullConfNode(bound, newNodes, leafNode, leafNode.getConfSearchNode());
                leafNode.markUpdated();
            }
            leafTime.stop();
            leafTimeAverage = leafTime.getTimeS();
            System.out.println("Processed " + numNodes + " leaves in " + leafTimeAverage + " seconds.");
            /*
             Experiment: assume we do enough work per conf to make additional parallelization unnecessary.
            if (bound.maxMinimizations < parallelism.numThreads/5)
                bound.maxMinimizations++;
             */
            bound.internalQueue.addAll(internalNodes);
        } else {
            numNodes = internalNodes.size();
            numInternals = numNodes;
            System.out.println("Processing " + numNodes + " internal nodes...");
            internalTime.reset();
            internalTime.start();
            for (MultiSequenceSHARKStarNode internalNode : internalNodes) {
                internalNode.checkDescendents(bound.sequence);
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
            internalTime.stop();
            internalTimeSum = internalTime.getTimeS();
            internalTimeAverage = internalTimeSum / Math.max(1, internalNodes.size());
            debugPrint("Internal node time :" + internalTimeSum + ", average " + internalTimeAverage);
            numInternalNodesProcessed += internalNodes.size();
            bound.leafQueue.addAll(leafNodes);
        }
        cleanupTime.start();
        loopCleanup(bound, newNodes, loopWatch, numNodes);
        cleanupTime.stop();

        if(writeTimes) {
            try {
                writer.write(String.format("%f, %f, %d, %f, %d, %f",
                        popQueuesTimer.getTimeS(),
                        internalTime.getTimeS(),
                        numInternals,
                        leafTime.getTimeS(),
                        numLeaves,
                        cleanupTime.getTimeS()
                ));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    protected void debugHeap(SHARKStarQueue queue) {
        queue.debugCheck(true);
    }


    boolean isStable(BigDecimal stabilityThreshold, SingleSequenceSHARKStarBound seqBound) {
        return seqBound.state.numEnergiedConfs <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(seqBound.state.getUpperBound(), stabilityThreshold);
    }


    protected void populateQueues(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> internalNodes,
                                  List<MultiSequenceSHARKStarNode> leafNodes, BigDecimal[] ZSums) {
        List<MultiSequenceSHARKStarNode> leftoverLeaves = new ArrayList<>();
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        int maxNodes = 1000;
        if (leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int) Math.floor(0.1 * leafTimeAverage / internalTimeAverage));
        while (!queue.isEmpty() && (bound.internalQueue.size() < maxNodes
                || (!bound.leafQueue.isEmpty() && MathTools.isGreaterThan(queue.peek().getErrorBound(bound.sequence),
                                                                            bound.leafQueue.peek().getErrorBound()))
                || (!bound.internalQueue.isEmpty() && MathTools.isGreaterThan(queue.peek().getErrorBound(bound.sequence),
                                                                            bound.internalQueue.peek().getErrorBound()))
                || bound.leafQueue.size() < bound.maxMinimizations)) {
            MultiSequenceSHARKStarNode curNode = queue.poll();
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
                String historyString = String.format("%s: previous lower bound %f, g score %f, hscore %f, f score %f corrected score %f, from %s",
                        node.confToString(), curNode.getConfLowerBound(bound.sequence), correctgscore, hscore, correctgscore+hscore, confCorrection, getStackTrace());
                curNode.setBoundsFromConfLowerAndUpperWithHistory(confCorrection, curNode.getConfUpperBound(bound.sequence), bound.sequence, historyString);
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
        ZSums[1] = fillListFromQueue(leafNodes, bound.leafQueue, bound.maxMinimizations, bound.sequence);
        if(!bound.internalQueue.isEmpty() &&
                MathTools.isLessThan(ZSums[0], bound.internalQueue.peek().getErrorBound(bound.sequence)))
            System.out.println("Should have used a node from the internal queue. How??");
        queue.addAll(leftoverLeaves);
    }

    private BigDecimal fillListFromQueue(List<MultiSequenceSHARKStarNode> list, Queue<MultiSequenceSHARKStarNode> queue, int max, Sequence seq) {
        BigDecimal sum = BigDecimal.ZERO;
        List<MultiSequenceSHARKStarNode> leftovers = new ArrayList<>();
        while (!queue.isEmpty() && list.size() < max ) {
            MultiSequenceSHARKStarNode curNode = queue.poll();
            if (correctedNode(leftovers, curNode, curNode.getConfSearchNode(), seq)) {
                BigDecimal diff = curNode.getUpperBound(seq).subtract(curNode.getLowerBound(seq));
                if(MathTools.isGreaterThan(diff, sum)) {
                    leftovers.remove(curNode);
                    queue.add(curNode);
                }
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
        for(int i =0; i<newNodes.size(); i++){
            MultiSequenceSHARKStarNode node = newNodes.get(i);
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
            node.setPartialConfLowerAndUpper(confCorrection, node.getPartialConfUpperBound());
            recordCorrection(oldg, confCorrection - oldg);
            String historyString = String.format("g score correction for %s: previous lower bound %f, g score %f, from %s",
                    node.confToString(), curNode.getConfLowerBound(seq), confCorrection, getStackTrace());
            curNode.setBoundsFromConfLowerAndUpperWithHistory(curNode.getConfLowerBound(seq) - oldg + confCorrection, curNode.getConfUpperBound(seq), seq, historyString);
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
                String historyString = "Error!";
                node.index(context.index);

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

                    historyString = String.format("%s: previous lower bound (none), g score %f, hscore %f, f score %f corrected score %f, from %s",
                            child.confToString(), diff, hdiff, diff+hdiff, confCorrection, getStackTrace());
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
                    historyString = String.format("%s: pairwise lower bound: %f, previous lower bound (none), correctedbound %f, from %s",
                            child.confToString(), lowerbound, confCorrection, getStackTrace());
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    confLowerBound = confCorrection;
                    child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                    confUpperBound = confRigid;
                    child.setPartialConfLowerAndUpper(confLowerBound, confUpperBound);
                    numConfsScored++;
                    progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
                }
                partialTime.stop();
                loopPartialTime += partialTime.getTimeS();

                child.index(context.index);
                SimpleConfSpace.Position designPos = confSpace.positions.get(nextPos);
                int nextDesignIndex = order.getNextPos(context.index, bound.seqRCs);
                SimpleConfSpace.Position nextDesignPos = null;
                if(nextDesignIndex >=0)
                    nextDesignPos = confSpace.positions.get(nextDesignIndex);
                MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild = curNode.makeOrUpdateChild(child, bound.sequence,
                        confLowerBound, confUpperBound, designPos, nextDesignPos);
                MultiSequenceSHARKStarNodeChild.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, bound.sequence, historyString);
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
        if (!subTreeRoot.hasChildren(precomputedSequence)){
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
        String historyString = "Error!!";

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
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();

            node.index(context.index);
            int nextPos = order.getNextPos(context.index, RCs);
            assert (!context.index.isDefined(nextPos));
            assert (context.index.isUndefined(nextPos));
            // score child nodes with tasks (possibly in parallel)
            List<MultiSequenceSHARKStarNode> children = new ArrayList<>();

            for (int nextRc : bound.seqRCs.get(nextPos)) {
                double resultingUpper = Double.MAX_VALUE;
                double resultingLower = -Double.MAX_VALUE;
                String historyString = "";

                node.index(context.index);
                Node child = node.assign(nextPos, nextRc);

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
                    //double lowerbound = minimizingEmat.confE(child.assignments);
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }
                    child.setBoundsFromConfLowerAndUpper(confLowerBound, confUpperbound);

                    if (debug) {
                        historyString = String.format("%s: previous lower bound (none), g score %f, hscore %f, f score %f corrected score %f, from %s",
                                node.confToString(), curNode.getConfLowerBound(bound.sequence), diff, hdiff, diff + hdiff, confCorrection, getStackTrace());
                    }
                    progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
                    resultingLower = confLowerBound;
                    resultingUpper= confUpperbound;
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
                    if (debug) {
                        historyString = String.format("%s: previous lower bound (none), confLower score %f, confCorrected score %f from %s",
                                node.confToString(), curNode.getConfLowerBound(bound.sequence), confLower, confCorrection, getStackTrace());
                    }
                    child.setBoundsFromConfLowerAndUpper(confCorrection, confRigid);
                    child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                    numConfsScored++;
                    progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
                    resultingLower= confCorrection;
                    resultingUpper= confRigid;
                }

                // Make the MSSHARK node

                child.index(context.index);
                SimpleConfSpace.Position designPos = confSpace.positions.get(nextPos);
                int nextDesignIndex = order.getNextPos(context.index, bound.seqRCs);
                SimpleConfSpace.Position nextDesignPos = null;
                if (nextDesignIndex >= 0)
                    nextDesignPos = confSpace.positions.get(nextDesignIndex);
                MultiSequenceSHARKStarNode newChild = curNode.makeOrUpdateChild(child,
                        bound.sequence, resultingLower, resultingUpper, designPos, nextDesignPos);
                newChild.setBoundsFromConfLowerAndUpperWithHistory(resultingLower,
                        resultingUpper, bound.sequence, historyString);
                //System.out.println("Created new child "+MultiSequenceSHARKStarNodeChild.toSeqString(bound.sequence));
                // collect the possible children
                if (newChild.getConfLowerBound(bound.sequence) < 0) {
                    children.add(newChild);

                }
                if (!newChild.isMinimized(bound.sequence)) {
                    newNodes.add(newChild);
                } else {
                    //newChild.computeEpsilonErrorBounds(bound.sequence);
                    //bound.addFinishedNode(newChild);
                    throw new RuntimeException();
                }

                //curNode.updateSubtreeBounds(bound.sequence);
                //printTree(bound.sequence, curNode);
                curNode.markUpdated();
            }
        }
    }

    public static boolean isDebugConf(int[] assignments) {
        return confMatch(debugConf, assignments);
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
        if(curNode.isMinimized(bound.sequence)) {
            bound.addFinishedNode(curNode);
            return;
        }
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
            node.setPartialConfLowerAndUpper(confCorrection, node.getPartialConfUpperBound());
            recordCorrection(oldg, confCorrection - oldg);
            String historyString = String.format("%s: correction from %f to %f, from ",
                    node.confToString(), oldg, confCorrection, getStackTrace());
            curNode.setBoundsFromConfLowerAndUpperWithHistory(confCorrection, curNode.getConfUpperBound(bound.sequence), bound.sequence, historyString);
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

                        //System.out.println("Minmized "+curNode.toSeqString(bound.sequence));
                        ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, curNode.getConfLowerBound(bound.sequence));
                        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                        Stopwatch correctionTimer = new Stopwatch().start();
                        if(batcher){
                            energyMatrixCorrector.scheduleEnergyCorrection(analysis, conf,
                                theBatcher);
                        }else{
                            computeEnergyCorrectionWithoutBatcher(analysis, conf, context.ecalc);
                        }

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

                        String historyString = String.format("minimimized %s to %s from %s",
                                node.confToString(), energy, getStackTrace());

                        curNode.setBoundsFromConfLowerAndUpperWithHistory(newConfLower, newConfUpper, bound.sequence, historyString);
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
                    seqBound.getLowerBound(), seqBound.getUpperBound(),
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
        public BatchCorrectionMinimizer batcher;
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

    private String getStackTrace() {
        return Arrays.stream(Thread.currentThread().getStackTrace())
                .map(Object::toString)
                .collect(Collectors.joining("\n"));
    }

    private void computeEnergyCorrectionWithoutBatcher(ConfAnalyzer.ConfAnalysis analysis, ConfSearch.ScoredConf conf,
                                         ConfEnergyCalculator ecalc) {
        if(conf.getAssignments().length < 3)
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
        for(int i = 0; i < sortedPairwiseTerms2.size(); i++)
        {
            TupE tupe = sortedPairwiseTerms2.get(i);
            double pairDiff = tupe.E;
            if(pairDiff < minDifference &&  maxDiff - pairDiff > threshhold)
                continue;
            maxDiff = Math.max(maxDiff, tupe.E);
            int pos1 = tupe.tup.pos.get(0);
            int pos2 = tupe.tup.pos.get(1);
            int localMinimizations = 0;
            for(int pos3 = 0; pos3 < diff.getNumPos(); pos3++) {
                if (pos3 == pos2 || pos3 == pos1)
                    continue;
                RCTuple tuple = makeTuple(conf, pos1, pos2, pos3);
                double tupleBounds = rigidEmat.getInternalEnergy(tuple) - minimizingEmat.getInternalEnergy(tuple);
                if(tupleBounds < triplethreshhold)
                    continue;
                minList.set(tuple.size()-1,minList.get(tuple.size()-1)+1);
                computeDifference(tuple, minimizingEcalc);
                localMinimizations++;
            }
            numPartialMinimizations+=localMinimizations;
            progress.reportPartialMinimization(localMinimizations, 0.0);
        }
        correctionTime.stop();
        //ecalc.tasks.waitForFinish();
    }

    private RCTuple makeTuple(ConfSearch.ScoredConf conf, int... positions) {
        RCTuple out = new RCTuple();
        for(int pos: positions)
            out = out.addRC(pos, conf.getAssignments()[pos]);
        return out;
    }
    private void computeDifference(RCTuple tuple, ConfEnergyCalculator ecalc) {
        computedCorrections = true;
        if(correctedTuples.contains(tuple.stringListing()))
            return;
        correctedTuples.add(tuple.stringListing());
        if(correctionMatrix.hasHigherOrderTermFor(tuple))
            return;
        minimizingEcalc.calcEnergyAsync(tuple, (minimizedTuple) -> {
            double tripleEnergy = minimizedTuple.energy;

            double lowerbound = minimizingEmat.getInternalEnergy(tuple);
            if (tripleEnergy - lowerbound > 0) {
                double correction = tripleEnergy - lowerbound;
                correctionMatrix.setHigherOrder(tuple, correction);
            }
            else
                System.err.println("Negative correction for "+tuple.stringListing());
        });
    }

    static class MultiSequenceState{
        long numEnergiedConfs = 0; // number of conformations fully minimized
        long numExpansions = 0; // number of internal nodes expanded
        long numPartialMinimizations; // number of partially minimized tuples

        double totalTimeEnergy = 0; // total time spent "energy-ing" conformations
        double totalTimeExpansion = 0; // total time spent expanding internal nodes
        double totalTimePartialMin = 0; // total time spent partially minimizing tuples

        long numRoundsEnergy = 0; // number of rounds of full minimization
        long numRoundsExpand = 0; // number of rounds of expansion
        long numRoundsPartialMin = 0; // number of rounds of partial minimization
    }

}
