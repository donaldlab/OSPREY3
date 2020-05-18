package edu.duke.cs.osprey.sharkstar_refactor;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseRigidGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.ScorerFactory;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.*;
import edu.duke.cs.osprey.energy.BatchCorrectionMinimizer;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.MARKStarProgress;
import edu.duke.cs.osprey.markstar.framework.StaticBiggestLowerboundDifferenceOrder;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.sharkstar.SHARKStarNodeScorer;
import edu.duke.cs.osprey.sharkstar.tools.SHARKStarEnsembleAnalyzer;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class MultiSequenceSHARKStarBound_refactor implements PartitionFunction {
    /**
     * A partition function calculator defined over a Multi-sequence conformation space
     *
     *
     * TODO: Implement updatePrecomputedNode
     * TODO: Work through tightenBoundinPhases
     *
     * TODO: When do I use loopTasks.waitForFinish?
     */

    // Debug variables
    public static final boolean debug = true;
    public static final boolean doExtraTupleCorrections = true;
    public static final boolean doUpperBoundCorrections = false;
    public final SHARKStarTreeDebugger pilotFish;
    public boolean profileOutput = false;

    // Settings
    protected double targetEpsilon = 1;
    private Status status = null;
    private MARKStarProgress progress;
    private String cachePattern = "NOT_INITIALIZED";
    private BigDecimal stabilityThreshold;


    // constructor variables
    private SimpleConfSpace confSpace;
    private EnergyMatrix rigidEmat;
    private EnergyMatrix minimizingEmat;
    public ConfEnergyCalculator minimizingEcalc;
    protected RCs fullRCs;
    protected Parallelism parallelism;

    // This will handle the parallelism?
    protected static TaskExecutor loopTasks;

    // Astar order to choose residue ordering
    public StaticBiggestLowerboundDifferenceOrder order;

    // things for the score context
    public ObjectPool<ScoreContext> contexts;
    private ScorerFactory gscorerFactory;
    private ScorerFactory rigidgscorerFactory;
    private ScorerFactory hscorerFactory;
    private ScorerFactory nhscorerFactory;

    // Matrix and corrector for energy corrections
    public final SimpleScoreCorrector lowerBoundCorrector;
    private final EnergyMatrixCorrector_refactor energyMatrixCorrector;

    // SHARK* specific variables
    public SHARKStarNode rootNode;
    private Sequence precomputedSequence;
    private MultiSequenceSHARKStarBound_refactor precomputedPfunc;
    public final BoltzmannCalculator bc;
    private double drillDownDifference; //The freeEnergy difference corresponding to target epsilon

    // Variables for upper bound corrections
    private TupleTrie<TupMapping> upperBoundCorrector;
    private RigidEmatFactory upperCorrectionEmatFactory;
    private EnergyMatrix upperCorrectionEmat;
    private SimpleConfSpace mutatedConfSpace;
    private boolean saveEPMOLsForMinimization = false;
    private Map<int[], EnergyCalculator.EnergiedParametricMolecule> epmolsMap;

    // Not sure where these should go
    private ConfAnalyzer confAnalyzer;
    private SHARKStarEnsembleAnalyzer ensembleAnalyzer;

    // Global tracking variables
    private double minimizationTimeTotal = 0.0;
    private double scoringTimeTotal = 0.0;
    private double correctionComputationTimeTotal = 0.0;

    private int numMinimizations = 0;
    private int numScores = 0;
    private int numCorrections = 0;

    private double leafTimeAverage;
    private double internalTimeAverage;
    private int numInternalNodesProcessed;

    public List<Integer> minList;

    // Corrections variables
    private boolean computedCorrections = false;
    private int numPartialMinimizations;
    private Set<String> correctedTuples;


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
    public MultiSequenceSHARKStarBound_refactor(SimpleConfSpace confSpace, EnergyMatrix rigidEmat, EnergyMatrix minimizingEmat,
                                                ConfEnergyCalculator minimizingConfEcalc, RCs rcs, Parallelism parallelism) {
        //, RigidEmatFactory upperCorrectionEmatFactory){

        // Set the basic needs
        this.confSpace = confSpace;
        this.rigidEmat = rigidEmat;
        this.minimizingEmat = minimizingEmat;
        this.minimizingEcalc = minimizingConfEcalc;
        this.fullRCs = rcs;

        // Initialize the correctionMatrix and the Corrector
        this.minList = new ArrayList<>(Collections.nCopies(rcs.getNumPos(), 0));
        this.lowerBoundCorrector = new SimpleScoreCorrector(confSpace.positions, MathTools.Optimizer.Maximize);
        //Intialize stuff for upperBoundCorrections
        //this.upperCorrectionEmatFactory = upperCorrectionEmatFactory;
        this.upperCorrectionEmatFactory = null;
        this.mutatedConfSpace = this.confSpace.makeCopy();
        this.upperBoundCorrector = new TupleTrie<>(confSpace.positions);

        // Set up the scoring machinery
        this.gscorerFactory = (emats) -> new PairwiseGScorer(emats);
        this.rigidgscorerFactory = (emats) -> new PairwiseRigidGScorer(emats);
        this.hscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, false);
        this.nhscorerFactory = (emats) -> new SHARKStarNodeScorer(emats, true);

        this.contexts = new ObjectPool<>((ignored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex(rcs.getNumPos());
            context.partialConfLBScorer = gscorerFactory.make(this.minimizingEmat);
            context.unassignedConfLBScorer = hscorerFactory.make(this.minimizingEmat);
            context.partialConfUBScorer = rigidgscorerFactory.make(this.rigidEmat);
            context.unassignedConfUBScorer = nhscorerFactory.make(this.rigidEmat); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;
            context.batcher = new BatchCorrectionMinimizer(minimizingConfEcalc, this.lowerBoundCorrector, this.minimizingEmat);

            return context;
        });

        // Boltzmann calculator for single sequence pfuncs to use
        bc = new BoltzmannCalculator(decimalPrecision);

        // Setting parallelism requires the ObjectPool being defined
        this.parallelism = parallelism;
        setParallelism(this.parallelism);

        // Initialize things for analyzing energies
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        ensembleAnalyzer = new SHARKStarEnsembleAnalyzer(minimizingEcalc, minimizingEmat);
        energyMatrixCorrector = new EnergyMatrixCorrector_refactor(this);

        progress = new MARKStarProgress(fullRCs.getNumPos());

        // Initialize debugger if necessary
        if (debug)
            pilotFish = new SHARKStarTreeDebugger(decimalPrecision);

        // No precomputed sequence means the "precomputed" sequence is empty
        this.precomputedSequence = confSpace.makeUnassignedSequence();
        pilotFish.setPrecomputedSequence(this.precomputedSequence);

    }

    public void setParallelism(Parallelism val) {

        if (val == null) {
            val = Parallelism.makeCpu(1);
        }

        parallelism = val;
        //loopTasks = minimizingEcalc.tasks;
        if (loopTasks == null)
            loopTasks = parallelism.makeTaskExecutor(null);
        contexts.allocate(parallelism.getParallelism());
    }

    @Override
    public void setReportProgress(boolean val) {

    }

    @Override
    public void setConfListener(ConfListener val) {

    }

    public void setCachePattern(String pattern) {
        this.cachePattern = pattern;
    }

    @Override
    public void setStabilityThreshold(BigDecimal threshold) {
        stabilityThreshold = threshold;
    }

    @Override
    public Status getStatus() {
        return null;
    }

    @Override
    public Values getValues() {
        return null;
    }

    @Override
    public int getParallelism() {
        return 0;
    }

    @Override
    public int getNumConfsEvaluated() {
        return 0;
    }

    public SimpleConfSpace getConfSpace() {
        return this.confSpace;
    }

    /**
     * Returns a wrapped pointer to this class, so that BBK* and MSK* can pretend they have single-sequence
     * partition functions.
     */
    public PartitionFunction getPartitionFunctionForSequence(Sequence seq) {
        SingleSequenceSHARKStarBound_refactor newBound = new SingleSequenceSHARKStarBound_refactor(this, seq, this.bc);
        newBound.init(null, null, targetEpsilon);
        System.out.println("Creating new pfunc for sequence " + seq);
        System.out.println("Full RCs: " + fullRCs);
        System.out.println("Sequence RCs: " + newBound.seqRCs);
        computeFringeForSequence(newBound, this.rootNode);
        newBound.state.setInitialBounds(
                newBound.calcEBound(n -> n.getFreeEnergyLB(newBound.sequence)),
                newBound.calcEBound(n -> n.getFreeEnergyUB(newBound.sequence))
        );
        if(debug){
            System.out.println("Fringe created, printing fringe level breakdown");
            System.out.println("Internal Queue:");
            SHARKStarQueueDebugger.printLevelBreakdown(newBound.internalQueue);
            System.out.println("Leaf Queue:");
            SHARKStarQueueDebugger.printLevelBreakdown(newBound.leafQueue);
        }
        // Wait for scoring to be done, if applicable
        loopTasks.waitForFinish();
        double boundEps = newBound.calcEpsilon();
        if (boundEps == 0) {
            System.err.println("Perfectly bounded sequence? how?");
        } else {
            System.out.println(String.format("Pfunc for %s created with epsilon of %.3f", seq, boundEps));
        }
        return newBound;
    }

    /**
     * init
     *
     * @param confSearch            We don't use this in SHARK*
     * @param numConfsBeforePruning We don't use this in SHARK*
     * @param targetEpsilon         The approximation error target
     */
    @Override
    public void init(ConfSearch confSearch, BigInteger numConfsBeforePruning, double targetEpsilon) {
        init(targetEpsilon);
    }

    /**
     * initialize partition function
     *
     * @param targetEpsilon Approximation error target
     */
    public void init(double targetEpsilon) {
        this.targetEpsilon = targetEpsilon;
        this.drillDownDifference = bc.freeEnergy(new BigDecimal(1 - targetEpsilon));
        this.status = Status.Estimating;
        makeAndScoreRootNode();
        if (precomputedPfunc == null)
            precomputeFlexible();
    }

    /**
     * initialize partition function
     * <p>
     * This should *only* be called if this is a flexible (non-mutable) partition function
     *
     * @param epsilon            Approximation error target
     * @param stabilityThreshold Minimum acceptable partition function value
     */
    private void initFlex(double epsilon, BigDecimal stabilityThreshold, boolean saveEPMOLsForMinimization) {
        this.targetEpsilon = epsilon;
        this.status = Status.Estimating;
        this.stabilityThreshold = stabilityThreshold;
        this.saveEPMOLsForMinimization = saveEPMOLsForMinimization;
        if (this.saveEPMOLsForMinimization)
            this.epmolsMap = new HashMap<>();
        makeAndScoreRootNode();
    }

    /**
     * Make and get starting bounds for the root node of the multi-sequence tree. Also initializes the residue ordering
     * so that using precomputed flexibility works properly.
     */
    public void makeAndScoreRootNode() {
        // Make the root node
        int[] rootConf = new int[confSpace.getNumPos()];
        Arrays.fill(rootConf, SHARKStarNode.Unassigned);
        this.rootNode = new SHARKStarNode(rootConf, 0, null);
        if (debug)
            pilotFish.setRootNode(this.rootNode);

        // Initialize residue ordering
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat));

        // Score the root node, being consistent even though we don't need to do it so fancy
        loopTasks.submit(
                () -> {
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        ScoreContext context = checkout.get();

                        ScoringResult result = new ScoringResult();
                        result.resultNode = rootNode;
                        rootNode.index(context.index);

                        // Force initialization of the residue ordering, kind of a hack,
                        // but necessary for flexible precomputation
                        order.getNextPos(context.index, fullRCs);

                        // Score the root
                        result.partialLB = context.partialConfLBScorer.calc(context.index, this.fullRCs);
                        result.partialUB = context.partialConfUBScorer.calc(context.index, this.fullRCs);
                        result.unassignLB = context.unassignedConfLBScorer.calc(context.index, this.fullRCs);
                        result.unassignUB = context.unassignedConfUBScorer.calc(context.index, this.fullRCs);
                        // Compute the node partition function error
                        result.score = bc.calc_lnZDiff(result.partialLB + result.unassignLB,
                                result.partialUB + result.unassignUB);

                        return result;
                    }
                },
                (result) -> {
                    if (!result.isValid())
                        throw new RuntimeException("Error in root node scoring");

                    result.resultNode.setPartialConfLB(result.partialLB);
                    result.resultNode.setPartialConfUB(result.partialUB);
                    result.resultNode.setUnassignedConfLB(result.unassignLB, this.precomputedSequence);
                    result.resultNode.setUnassignedConfUB(result.unassignUB, this.precomputedSequence);
                    result.resultNode.setScore(result.score, this.precomputedSequence);
                }
        );
        loopTasks.waitForFinish();

        System.out.println(String.format("Initial root free energy bounds: [%.3f, %.3f]",
                this.rootNode.getFreeEnergyLB(this.precomputedSequence),
                this.rootNode.getFreeEnergyUB(this.precomputedSequence)));
    }

    /**
     * Precompute the partition function for the flexible residues
     * <p>
     * We are recomputing the energy matrices here because there's some
     * stupidly nontrivial parallel array mapping to be done to ensure that
     * the extensive energy calculation features we rely on are
     * computing the right energy for the flexible conf space.
     * <p>
     * I'm sure this could be optimized, but I doubt that it's worth it,
     * since this only happens once per state.
     */
    public void precomputeFlexible() {
        // Make a copy of the confSpace without mutable residues
        SimpleConfSpace flexConfSpace = confSpace.makeFlexibleCopy();
        Sequence unassignedFlex = flexConfSpace.makeUnassignedSequence();

        // Sometimes our designs don't have immutable residues on one side.
        if (flexConfSpace.positions.size() > 0) {
            System.out.println("Making flexible confspace bound...");

            // Make an all-new SHARKStarBound for the flexible partition function
            RCs flexRCs = new RCs(flexConfSpace);
            ConfEnergyCalculator flexMinimizingConfECalc = FlexEmatMaker.makeMinimizeConfEcalc(flexConfSpace,
                    this.parallelism);
            ConfEnergyCalculator rigidConfECalc = FlexEmatMaker.makeRigidConfEcalc(flexMinimizingConfECalc);
            EnergyMatrix flexMinimizingEmat = FlexEmatMaker.makeEmat(flexMinimizingConfECalc, "minimized", cachePattern + ".flex");
            EnergyMatrix flexRigidEmat = FlexEmatMaker.makeEmat(rigidConfECalc, "rigid", cachePattern + ".flex");

            // Construct the bound
            MultiSequenceSHARKStarBound_refactor precompFlex = new MultiSequenceSHARKStarBound_refactor(
                    flexConfSpace, flexRigidEmat, flexMinimizingEmat,
                    flexMinimizingConfECalc, flexRCs, this.parallelism);
            // Set settings
            precompFlex.setCachePattern(cachePattern);
            // initialize
            precompFlex.initFlex(this.targetEpsilon, this.stabilityThreshold, MultiSequenceSHARKStarBound_refactor.doUpperBoundCorrections);
            PartitionFunction flexBound = precompFlex.getPartitionFunctionForSequence(unassignedFlex);
            flexBound.compute();
            precompFlex.printEnsembleAnalysis();
            // Copy over the minimized fragments
            this.epmolsMap = precompFlex.epmolsMap;
            processPrecomputedFlex(precompFlex);
            if(doUpperBoundCorrections && !saveEPMOLsForMinimization){
                this.upperCorrectionEmat = this.upperCorrectionEmatFactory.make(mutatedConfSpace);
                try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                    ScoreContext context = checkout.get();
                    context.upperCorrectionGScorer = this.rigidgscorerFactory.make(this.upperCorrectionEmat);
                    context.upperCorrectionHScorer = this.nhscorerFactory.make(this.upperCorrectionEmat);
                }
            }
            System.out.println(String.format("Precomputation of flexible residues resulted in %d partial minimizations", lowerBoundCorrector.getNumCorrections()));
        }
    }

    /**
     * Cannibalizes a precomputed pfunc to start this pfunc
     *
     * @param precomputedFlex The precomputed partition function
     */
    private void processPrecomputedFlex(MultiSequenceSHARKStarBound_refactor precomputedFlex) {
        precomputedPfunc = precomputedFlex;
        //precomputedRootNode = precomputedFlex.rootNode;
        this.precomputedSequence = precomputedFlex.confSpace.makeWildTypeSequence();
        updatePrecomputedConfTree(precomputedFlex.rootNode);
        mergeCorrections(precomputedFlex.lowerBoundCorrector, genConfSpaceMapping());

        // Fix order issues. Hacky hack
        ConfIndex rootIndex = new ConfIndex(fullRCs.getNumPos());
        this.rootNode.index(rootIndex);
        this.order.updateForPrecomputedOrder(precomputedFlex.order, rootIndex, this.fullRCs, genConfSpaceMapping());
    }

    /**
     * Makes the precomputed confTree consistent with the full confSpace
     * <p>
     * When we precompute flexible residues, we will have a tree that is for a flexible confspace.
     * However, when we want to compute for mutable residues, we need to extend the length of assignments in our tree
     *
     * @param precomputedRootNode The root of the precomputed partition function
     *                            <p>
     *                            TODO: Make this work even if we stop storing root nodes
     */
    private void updatePrecomputedConfTree(SHARKStarNode precomputedRootNode) {
        int[] permutationArray = genConfSpaceMapping();
        updatePrecomputedNode(precomputedRootNode, permutationArray, this.confSpace.getNumPos());
        this.rootNode = precomputedRootNode;
    }

    /**
     * Recursive helper function for updatePrecomputedConfTree
     *
     * @param node        The current node to update
     * @param permutation The permutation matrix for mapping the precomputed RCs to the current RCs
     * @param size        The size of the new confSpace
     */
    private void updatePrecomputedNode(SHARKStarNode node, int[] permutation, int size) {
        // Store old assignments for epmolsMap
        int[] oldAssignments = node.getAssignments();
        // permute the assignments to the new confspace
        node.makeNodeCompatibleWithConfSpace(permutation, size);

        //set corrections if node is minimized
        if (node.isMinimized()) {
            RCTuple nodeTup = new RCTuple(node.getAssignments());
            lowerBoundCorrector.insertCorrection(new TupE(nodeTup, node.getMinE()
                    - node.getPartialConfLB()));
            if(doUpperBoundCorrections && !saveEPMOLsForMinimization && epmolsMap.containsKey(oldAssignments)){
                // Add minimized fragment to mutated confspace
                RCTuple mappedTup = mutatedConfSpace.addResidueConfsFromEPMOL(nodeTup, epmolsMap.get(oldAssignments));
                // Insert the correction
                this.upperBoundCorrector.insert(new TupMapping(nodeTup, mappedTup));
            }
        }
        /*
        Reset minimized energies and unassignedConf bounds, since these will no longer be valid.
        Partial conf bounds should still be valid
        Also reset upperbound corrections, which will be invalid
         */
        node.setIsMinimized(false);
        node.setMinE(Double.NaN);
        node.setUnassignedConfLB(Double.NaN, precomputedSequence);
        node.setUnassignedConfUB(Double.NaN, precomputedSequence);
        node.unsetHOTCorrectionUB();

        // Recurse
        if (!node.getChildren().isEmpty()) {
            for (SHARKStarNode child : node.getChildren()) {
                updatePrecomputedNode(child, permutation, size);
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

    /**
     * Takes partial minimizations from the precomputed correctionMatrix, maps them to the new confspace, and
     * stores them in this correctionMatrix
     */
    public void mergeCorrections(SimpleScoreCorrector precomputedCorrections, int[] confSpacePermutation) {
        List<TupE> corrections = precomputedCorrections.getAllCorrections().stream()
                .map((tup) -> tup.permute(confSpacePermutation))
                .collect(Collectors.toList());
        if (corrections.size() != 0) {
            int TestNumCorrections = corrections.size();
            this.lowerBoundCorrector.insertAllCorrections(corrections);
        } else
            System.out.println("No corrections to insert");
    }

    /**
     * Gets children that are compatible with a given sequence rCs object
     *
     * @param node      The node for which to get children
     * @param seqRCs    The RCs for each sequence position
     * @return
     *
     * TODO: Determine whether it is better to have a datastructure in the node that keeps track of this info
     */
    public List<SHARKStarNode> getChildrenCompatibleWithSeqRCs(SHARKStarNode node, RCs seqRCs) {
        ArrayList<SHARKStarNode> compatibleChildren = new ArrayList<>();
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            node.index(context.index);
            // It shouldn't matter whether we use fullRCs or seqRCs here, since this position should already be defined
            int nextPos = order.getNextPos(context.index, fullRCs);
            // If the next pos doesn't exist, just return an empty list
            if(nextPos == -1){
                return compatibleChildren;
            }
            // Get the RCs allowed at the given Pos
            int[] allowedRCsAtPos = seqRCs.get(nextPos);
            // Check children using loops, since it's probably fastest
            for (SHARKStarNode child : node.getChildren()) {
                int childRCAtPos = child.getAssignments()[nextPos];
                for (int rc : allowedRCsAtPos) {
                    if (childRCAtPos == rc) {
                        compatibleChildren.add(child);
                    }
                }
            }
            return compatibleChildren;
        }
    }

    /**
     * Compute the portion of the multi-sequence tree fringe nodes that are compatible with a given sequence
     * @param bound The single sequence bound
     * @param node  The current node, used for recursion
     *
     * TODO: Make this just look through the fringe instead of starting from the root
     */
    public void computeFringeForSequence(SingleSequenceSHARKStarBound_refactor bound, SHARKStarNode node) {
        RCs rcs = bound.seqRCs;
        // Get a list of sequence-compatible children for this node
        List<SHARKStarNode> compatibleChildren = getChildrenCompatibleWithSeqRCs(node, rcs);
        // If no compatible children, this could be a fringe node
        if (compatibleChildren.isEmpty()) {
            // Check to see if there is currently an unassigned energy calculated for this sequence
            if (!node.getUnassignedConfLB().containsKey(bound.sequence) &&
                    !node.getUnassignedConfUB().containsKey(bound.sequence)) {
                // If there is not, score the node
                scoreNodeForSeq(node, bound.sequence, bound.seqRCs);
                //TODO: determine whether I'll have thread issues by not waiting till this is finished to add the node to fringe
            }
            loopTasks.waitForFinish();
            // add to sequence fringe nodes
            if (node.getLevel() < bound.seqRCs.getNumPos()) {
                bound.internalQueue.put(node);
            } else {
                bound.leafQueue.put(node);
            }
        } else {
            // If there are compatible children, recurse
            for (SHARKStarNode child : compatibleChildren) {
                computeFringeForSequence(bound, child);
            }
        }

    }

    /**
     * Use the Score context to score a single node (possibly parallel)
     *
     * @param node   The node to score
     * @param seqRCs The seqRCs to score over
     */
    public void scoreNodeForSeq(SHARKStarNode node, Sequence seq, RCs seqRCs) {
        loopTasks.submit(
                () -> {
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        // Timing
                        /* Note: this is probably inefficient to create new stopwatches, but they aren't threadsafe,
                        so I'm doing this for now
                         */
                        Stopwatch scoreWatch = new Stopwatch();
                        scoreWatch.start();

                        ScoreContext context = checkout.get();

                        ScoringResult result = new ScoringResult();
                        result.resultNode = node;
                        node.index(context.index);

                        // Score the node
                        //TODO: move whatever possible back to calcDifferential
                        result.partialLB = context.partialConfLBScorer.calc(context.index, seqRCs);
                        result.partialUB = context.partialConfUBScorer.calc(context.index, seqRCs);
                        result.unassignLB = context.unassignedConfLBScorer.calc(context.index, seqRCs);
                        result.unassignUB = context.unassignedConfUBScorer.calc(context.index, seqRCs);
                        // Compute the node partition function error
                        result.score = bc.calc_lnZDiff(result.partialLB + result.unassignLB,
                                result.partialUB + result.unassignUB);

                        scoreWatch.stop();
                        result.time = scoreWatch.getTimeS();

                        return result;
                    }
                },
                (result) -> {
                    if (!result.isValid())
                        throw new RuntimeException(String.format("Error in node scoring for %s",
                                this.confSpace.formatConf(result.resultNode.getAssignments())));

                    result.resultNode.setPartialConfLB(result.partialLB);
                    result.resultNode.setPartialConfUB(result.partialUB);
                    result.resultNode.setUnassignedConfLB(result.unassignLB, seq);
                    result.resultNode.setUnassignedConfUB(result.unassignUB, seq);
                    result.resultNode.setScore(result.score, seq);

                    synchronized (this) {
                        scoringTimeTotal += result.time;
                        numScores += 1;
                    }
                }
        );
    }

    /**
     * Make a child node and score it
     *
     * @param parent  The parent of this child
     * @param nextPos The residue position to assign for the child
     * @param nextRC  The residue conformation to assign to the child
     * @param seq     The sequence for the child
     * @param seqRCs  The seqRCs to score over
     */
    public void makeAndScoreChildForSeq(SHARKStarNode parent, int nextPos, int nextRC, Sequence seq, RCs seqRCs, List<SHARKStarNode> children) {

        loopTasks.submit(
                () -> {
                    ScoringResult result = new ScoringResult();
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        // Timing
                        /* Note: this is probably inefficient to create new stopwatches, but they aren't threadsafe,
                        so I'm doing this for now
                         */
                        Stopwatch scoreWatch = new Stopwatch();
                        scoreWatch.start();

                        ScoreContext context = checkout.get();

                        parent.index(context.index);
                        result.resultNode = parent.assign(nextPos, nextRC);

                        // Score the child node
                        //TODO move back to calcDifferential if possible
                        result.partialLB = context.partialConfLBScorer.calcDifferential(context.index, seqRCs, nextPos, nextRC);
                        result.partialUB = context.partialConfUBScorer.calcDifferential(context.index, seqRCs, nextPos, nextRC);
                        result.resultNode.index(context.index); // We don't implement calcDifferential for the SHARKStarNodeScorer yet, so probably faster to do this
                        result.unassignLB = context.unassignedConfLBScorer.calc(context.index, seqRCs);
                        result.unassignUB = context.unassignedConfUBScorer.calc(context.index, seqRCs);

                        scoreWatch.stop();
                        result.time = scoreWatch.getTimeS();

                        // Try to apply partial minimization correction to child (Lowerbound)
                        // If the parent correction is larger, just use the parent correction, since it must be valid for lowerbound
                        result.HOTCorrectionLB = Math.max(lowerBoundCorrector.getCorrection(context.index), parent.getHOTCorrectionLB());

                        if(doUpperBoundCorrections && !saveEPMOLsForMinimization && upperBoundCorrector.getEntries(new RCTuple(context.index.makeConf())).size()>0) {
                            // Try to apply correction to child (upperbound)
                            RCTuple mappedTup = upperBoundCorrector.getEntries(new RCTuple(context.index.makeConf())).get(0).mapTup;
                            mappedTup.pasteToIndex(context.index);
                            double correctedPartialUB = context.upperCorrectionGScorer.calc(context.index, seqRCs);
                            double correctedUnnasignedUB = context.upperCorrectionHScorer.calc(context.index, seqRCs);
                            double upperCorrectionSize = correctedPartialUB + correctedUnnasignedUB - (result.partialUB + result.unassignUB);
                            // Filter out positive upper bound corrections
                            result.HOTCorrectionUB = Math.min(0, upperCorrectionSize);
                        }else{
                            result.HOTCorrectionUB=0;
                        }

                        // Compute the node partition function error
                        result.score = bc.calc_lnZDiff(result.partialLB + result.unassignLB + result.HOTCorrectionLB,
                                result.partialUB + result.unassignUB + result.HOTCorrectionUB);

                        }
                    return result;
                },
                (result) -> {
                    if (!result.isValid())
                        throw new RuntimeException(String.format("Error in node scoring for %s",
                                this.confSpace.formatConf(result.resultNode.getAssignments())));

                        result.resultNode.setPartialConfLB(result.partialLB);
                        result.resultNode.setPartialConfUB(result.partialUB);
                        result.resultNode.setUnassignedConfLB(result.unassignLB, seq);
                        result.resultNode.setUnassignedConfUB(result.unassignUB, seq);
                        result.resultNode.setHOTCorrectionLB(result.HOTCorrectionLB);
                        if (doUpperBoundCorrections && result.HOTCorrectionUB < 0)
                            result.resultNode.setHOTCorrectionUB(result.HOTCorrectionUB);
                        result.resultNode.setScore(result.score, seq);

                    synchronized (this) {
                        scoringTimeTotal += result.time;
                        numScores += 1;
                        children.add(result.resultNode);
                    }

                }
        );
    }

    /**
     * Use the scoreContext to minimize a single node with tasks
     *
     * @param node The node to minimize
     * @param seq  The sequence whose score we will compare to the minimized energy
     *             <p>
     *             TODO: remove the dependency on sequence, it's not necessary
     */
    public void minimizeNodeForSeq(SHARKStarNode node, Sequence seq) {
        loopTasks.submit(
                () -> {
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        //Timing
                        Stopwatch minimizationTimer = new Stopwatch();
                        Stopwatch correctionTimer = new Stopwatch();
                        minimizationTimer.start();

                        ScoreContext context = checkout.get();
                        MinimizationResult result = new MinimizationResult();
                        result.resultNode = node;
                        node.index(context.index);

                        ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.getAssignments(), node.getFreeEnergyLB(seq));
                        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
                        node.setMinE(analysis.epmol.energy);
                        result.minimizedEnergy = analysis.epmol.energy;
                        node.setScore(0.0, seq);
                        if(this.saveEPMOLsForMinimization)
                            result.epmol = analysis.epmol;

                        minimizationTimer.stop();
                        if (doExtraTupleCorrections) {
                            correctionTimer.start();
                            energyMatrixCorrector.computeEnergyCorrection(analysis, conf, 1.0,
                                    context.batcher);
                            correctionTimer.stop();
                        }
                        result.minimizationTime = minimizationTimer.getTimeS();
                        result.correctionTime = correctionTimer.getTimeS();
                    /*
                    //TODO: add corrections back in
                    correctionTimer.start();

                    energyMatrixCorrector.computeEnergyCorrection(analysis, conf, bound.getSequenceEpsilon(),
                            context.batcher);

                     */

                        return result;
                    }
                },
                (result) -> {
                    if (!result.isValid())
                        throw new RuntimeException(String.format("Error in node scoring for %s",
                                this.confSpace.formatConf(result.resultNode.getAssignments())));

                    System.out.println(String.format("Minimized %s --> %.3f in [%.3f, %.3f]",
                            result.resultNode.confToString(),
                            result.minimizedEnergy,
                            result.resultNode.getFreeEnergyLB(seq),
                            result.resultNode.getFreeEnergyUB(seq)
                    ));

                    result.resultNode.setMinE(result.minimizedEnergy); // set the energy
                    result.resultNode.setIsMinimized(true); // set the minimized flag
                    result.resultNode.setScore(0.0, seq); // Since the node is minimized, the score (error) is 0

                    minList.set(result.resultNode.getAssignments().length - 1, minList.get(result.resultNode.getAssignments().length - 1) + 1);


                    synchronized (this) {
                        minimizationTimeTotal += result.minimizationTime;
                        correctionComputationTimeTotal += result.correctionTime;
                        numMinimizations += 1;
                        // Save epmol if we are keeping it
                        if(result.epmol != null){
                            epmolsMap.put(result.resultNode.getAssignments(), result.epmol);
                        }
                    }
                });
    }

    @Override
    public void compute(int maxNumConfs) {
        throw new UnsupportedOperationException("Do not try to run Multisequence SHARK* bounds directly. Call " +
                "getPartitionFunctionForSequence() and use the generated single sequence bound.");
    }

    /**
     * Compute the partition function for a given sequence using the multi-sequence tree. Halts once it reaches
     * maxNumConfs or when it reaches this.targetEpsilon, whichever comes first.
     *
     * @param maxNumConfs The maximum number of conformations to evaluate (minimize?)
     * @param bound       The SingleSequence partition function
     */
    public void computeForSequence(int maxNumConfs, SingleSequenceSHARKStarBound_refactor bound) {
        if (debug){
            pilotFish.setRootNode(this.rootNode);
        }
        Stopwatch computeWatch = new Stopwatch().start();

        System.out.println("Tightening bound for " + bound.sequence);
        debugPrint("Num conformations: " + bound.numConformations);

        double lastEps = bound.calcEpsilon();
        if (lastEps == 0)
            System.err.println("ERROR: Computing for a partition function with epsilon=0");

        /*
        //TODO: Implement this optimization
        if (!bound.nonZeroLower()) {
            runUntilNonZero(sequenceBound);
            sequenceBound.updateBound();
        }
         */

        //TODO: make this pay attention to maxNumConfs by getting a workDone() method
        // and make this pay attention to stability calculations
        while (lastEps > targetEpsilon
            //workDone() - previousConfCount < maxNumConfs &&
            //isStable(stabilityThreshold, sequenceBound.sequence
        ) {
            debugPrint("Tightening from epsilon of " + lastEps);
            if (debug) {
                //pilotFish.travelTree(bound.sequence);
            }
            // run method
            //tightenBoundInPhases(bound);
            //simpleTightenBound(bound);
            try {
                moreComplicatedTightenBound(bound);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

            // do some debug checks
            double newEps = bound.calcEpsilon();
            debugPrint("Errorbound is now " + newEps);
            debugPrint("Bound reduction: " + (lastEps - newEps));
            if (lastEps < newEps
                //|| bound.errors()
            ) {
                System.err.println("Error. Bounds got looser.");
                pilotFish.travelTree(bound.sequence);
                System.exit(-1);
            }
            lastEps = newEps;
        }
        //TODO: make stability threshold stuff
        /*
        if (!isStable(stabilityThreshold, sequenceBound.sequence))
            sequenceBound.setStatus(Status.Unstable);

         */
        loopTasks.waitForFinish();
        minimizingEcalc.tasks.waitForFinish();
        computeWatch.stop();

        System.out.println(String.format("Finished computation in %.3f seconds", computeWatch.getTimeS()));
        System.out.println(String.format("Final E bounds [%.3f, %.3f]]",
                bound.getELB(),
                bound.getEUB()
        ));
        System.out.println(String.format("Final Z bounds [%9.10e, %9.10e]]",
                bound.getZLB(),
                bound.getZUB()
        ));
        System.out.println(String.format("Direct computation of Z bounds [%9.10e, %9.10e], eps %.9f",
                bound.calcZLBDirect().doubleValue(),
                bound.calcZUBDirect().doubleValue(),
                bound.calcEpsilonDirect()));

        if (SHARKStarQueueDebugger.hasDuplicateNodes(bound.internalQueue)) {
            System.out.println("Internal Queue has duplicates");
        } else {
            System.out.println("Internal Queue does not have duplicates");
        }

        SHARKStarQueueDebugger.printLevelBreakdown(bound.internalQueue);

        System.out.println(String.format("Minimized %d conformations.",numMinimizations));

        if(debug) {
            pilotFish.travelTree(bound.sequence);
            //pilotFish.printLeaves(bound.sequence);
        }

        System.out.println(String.format("For internal queue:\n\t Z reduction from corrections: %1.5e, [lower, upper] = [%1.5e, %1.5e]"+
                "\n\t Z reduction from minimizations: %1.5e",
                SHARKStarQueueDebugger.getZErrorReductionFromCorrections(bound.internalQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromLowerCorrections(bound.internalQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromUpperCorrections(bound.internalQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromMinimizations(bound.internalQueue, bound.sequence).doubleValue()
                ));
        System.out.println(String.format("For leaf queue:\n\t Z reduction from corrections: %1.5e, [lower, upper] = [%1.5e, %1.5e]"+
                        "\n\t Z reduction from minimizations: %1.5e",
                SHARKStarQueueDebugger.getZErrorReductionFromCorrections(bound.leafQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromLowerCorrections(bound.leafQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromUpperCorrections(bound.leafQueue, bound.sequence).doubleValue(),
                SHARKStarQueueDebugger.getZErrorReductionFromMinimizations(bound.leafQueue, bound.sequence).doubleValue()
        ));
    }

    /**
     * Uses the simplest possible scheme to tighten the bounds on the partition function
     *
     * @param bound A single-sequence partition function
     */
    private void moreComplicatedTightenBound(SingleSequenceSHARKStarBound_refactor bound) throws InterruptedException {
        int numConfsToProcess = 10;
        AtomicInteger numConfsProcessed = new AtomicInteger();
        boolean leavesExist = false;
        double bestInternalScore = Double.NEGATIVE_INFINITY;
        double bestLeafScore = Double.NEGATIVE_INFINITY;


        System.out.println(String.format("Internal Queue: %d nodes, Leaf Queue: %d nodes",bound.internalQueue.size(), bound.leafQueue.size()));

        while (numConfsProcessed.get() < numConfsToProcess) {
            SHARKStarNode node;
            List<SHARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());

            Step step = Step.None;
            synchronized(this) {
                leavesExist = !bound.leafQueue.isEmpty();
                if (leavesExist)
                    bestLeafScore = bound.leafQueue.peek().getScore(bound.sequence);
                else
                    bestLeafScore = Double.NEGATIVE_INFINITY;
            }
            synchronized(this) {
                if (!bound.internalQueue.isEmpty())
                    bestInternalScore = bound.internalQueue.peek().getScore(bound.sequence);
                else
                    bestInternalScore = Double.NEGATIVE_INFINITY;
            }
            synchronized(this){
                System.out.println(String.format("Leaf queue score %.3f contains %d nodes, internal queue score %.3f contains %d nodes",
                        bestLeafScore, bound.leafQueue.size(), bestInternalScore, bound.internalQueue.size()));

                if (!leavesExist || bestInternalScore > bestLeafScore) {
                    step = Step.Score;
                }else {
                    step = Step.Energy;
                }
            }

            switch(step) {
                case Score: {
                    double finalBestLeafScore = bestLeafScore;
                    loopTasks.submit( () -> processPartialConfNodes(bound, 200, finalBestLeafScore),

                            (result) -> {
                                synchronized(this){
                                    numConfsProcessed.addAndGet(result.numNodesCreated);
                                }
                                synchronized (bound){
                                    if(result.numNodesCreated > 0){
                                        bound.state.addDiffToLB(bc.logSumExp(result.lowerBoundDeltas));
                                        bound.state.subtractDiffFromUB(bc.logSumExp(result.upperBoundDeltas));
                                    }
                                }
                            }
                    );
                    break;
                }
                case Energy: {
                    loopTasks.submit( () -> minimizeNode(bound),

                            (result) -> {
                                synchronized (this) {
                                    //minList.set(result.resultNode.getAssignments().length - 1, minList.get(result.resultNode.getAssignments().length - 1) + 1);
                                    minimizationTimeTotal += result.minimizationTime;
                                    correctionComputationTimeTotal += result.correctionTime;
                                    numMinimizations += 1;
                                    numConfsProcessed.addAndGet(1);

                                }
                                synchronized (bound){
                                    bound.state.addDiffToLB(result.lowerBoundDelta);
                                    bound.state.subtractDiffFromUB(result.upperBoundDelta);
                                }
                            }
                    );
                    break;
                }
                case None: {
                    break;
                }
            }
        }
        loopTasks.waitForFinish();
    }

    /**
     * Uses the simplest possible scheme to tighten the bounds on the partition function
     *
     * @param bound A single-sequence partition function
     */
    private void simpleTightenBound(SingleSequenceSHARKStarBound_refactor bound) {
        int numConfsToProcess = 1000;
        AtomicInteger numConfsProcessed = new AtomicInteger();
        boolean leavesExist = false;
        double bestInternalScore = Double.NEGATIVE_INFINITY;
        double bestLeafScore = Double.NEGATIVE_INFINITY;

        System.out.println(String.format("Internal Queue: %d nodes, Leaf Queue: %d nodes",bound.internalQueue.size(), bound.leafQueue.size()));

        while (numConfsProcessed.get() < numConfsToProcess) {
            SHARKStarNode node;
            List<SHARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>());
            leavesExist = !bound.leafQueue.isEmpty();
            bestInternalScore = bound.internalQueue.peek().getScore(bound.sequence);
            if(leavesExist)
                bestLeafScore = bound.leafQueue.peek().getScore(bound.sequence);

            if (!leavesExist || bestInternalScore > bestLeafScore) {
                node = bound.internalQueue.poll();
                processPartialConfNode(bound, newNodes, node);
            } else {
                node = bound.leafQueue.poll();
                processFullConfNode(bound, newNodes, node);
            }

            loopTasks.waitForFinish();
            for (SHARKStarNode newNode : newNodes) {
                if (newNode.getLevel() < bound.seqRCs.getNumPos()) {
                    bound.internalQueue.put(newNode);
                } else {
                    bound.leafQueue.put(newNode);
                }
            }

            numConfsProcessed.getAndIncrement();
        }
    }

    /**
     * Tighten bounds on a pfunc by following and choosing from the following steps:
     * <p>
     * 1.   Determine which nodes we want to process
     * 2.   Choose from choices A or B below
     * A) Process leaf nodes
     * B) Process internal nodes
     * 3.   Cleanup
     *
     * @param bound A single-sequence partition function
     */
    private void tightenBoundInPhases(SingleSequenceSHARKStarBound_refactor bound) {
        assert (!bound.isEmpty());

        // Initialize lists to track nodes
        List<SHARKStarNode> internalNodes = new ArrayList<>(); // List of internal nodes to score and correct
        List<SHARKStarNode> leafNodes = new ArrayList<>();     // List of leaf nodes to correct or minimize
        List<SHARKStarNode> newNodes = Collections.synchronizedList(new ArrayList<>()); // List of nodes to put back into the single-sequence bound fringe queue

        // Initialize tracking variables
        // NB: each of these variables tracks the Z *error* not the Z contribution
        BigDecimal internalZ = BigDecimal.ONE;
        BigDecimal leafZ = BigDecimal.ONE;
        BigDecimal[] ZSums = new BigDecimal[]{internalZ, leafZ};
        int numNodes = 0;

        // Initialize stopwatches and variables for timing
        Stopwatch loopWatch = new Stopwatch();
        loopWatch.start();
        Stopwatch internalTime = new Stopwatch();
        Stopwatch leafTime = new Stopwatch();
        double leafTimeSum = 0;
        double internalTimeSum = 0;

        // Record the current information before taking nodes out
        double startingELB = bound.calcEBound(e -> e.getFreeEnergyLB(bound.sequence));
        double startingEUB = bound.calcEBound(e -> e.getFreeEnergyUB(bound.sequence));
        double startingEps = bound.calcEpsilon();

        System.out.println(String.format("Current overall error bound: %.10f, spread of [%.3f, %.3f]",
                startingEps,
                startingELB,
                startingEUB));
        /*
        Read the bound fringeNodes and decide which nodes we want to process
         */
        populateQueues(bound, internalNodes, leafNodes, ZSums);

        // debugging the populateQueues method
        internalZ = ZSums[0];
        leafZ = ZSums[1];
        // If we don't have any nodes to process, try again
        if (leafNodes.isEmpty() && internalNodes.isEmpty()) {
            System.out.println("Nothing was populated?");
            populateQueues(bound, internalNodes, leafNodes, ZSums);
        }
        if (MathTools.isRelativelySame(internalZ, leafZ, PartitionFunction.decimalPrecision, 1e-3)
                && MathTools.isRelativelySame(leafZ, BigDecimal.ZERO, PartitionFunction.decimalPrecision, 1e-3)) {
            //pilotFish.travelTree(bound.sequence);
            //printTree(bound.sequence, rootNode);
            System.out.println("This is a bad time.");
            populateQueues(bound, new ArrayList<>(), new ArrayList<>(), ZSums);
        }
        System.out.println(String.format("Z Comparison: %12.6e, %12.6e", internalZ, leafZ));
        if (!bound.internalQueue.isEmpty() &&
                MathTools.isLessThan(internalZ, bc.calc(bound.internalQueue.peek().getFreeEnergyLB(bound.sequence))))
            System.out.println("Should have used a node from the internal queue. How??");

        /*
            If the internal nodes have less error than the leaf nodes, tighten bounds on leaf nodes
         */
        if (MathTools.isLessThan(internalZ, leafZ)) {
            numNodes = leafNodes.size();
            System.out.println("Processing " + numNodes + " leaf nodes...");
            leafTime.reset();
            leafTime.start();
            for (SHARKStarNode leafNode : leafNodes) {
                //debugPrint("Processing Node: " + leafNode.toSeqString(bound.sequence));
                processFullConfNode(bound, newNodes, leafNode);
                // do we need a markupdated method?
            }
            loopTasks.waitForFinish();
            leafTime.stop();
            leafTimeAverage = leafTime.getTimeS();
            System.out.println("Processed " + numNodes + " leaves in " + leafTimeAverage + " seconds.");
            /*
             Experiment: assume we do enough work per conf to make additional parallelization unnecessary.
            if (bound.maxMinimizations < parallelism.numThreads/5)
                bound.maxMinimizations++;
             */
            bound.internalQueue.addAll(internalNodes);
        /*
            If the leaf nodes have less error than the internal nodes, tighten bounds on internal nodes
         */
        } else {
            numNodes = internalNodes.size();
            System.out.println("Processing " + numNodes + " internal nodes...");
            internalTime.reset();
            internalTime.start();
            for (SHARKStarNode internalNode : internalNodes) {
                //TODO: test to make sure there are no children compatible

                System.out.println(String.format("Internal node: %s with free energy bounds [%.3f, %.3f]",
                        Arrays.toString(internalNode.getAssignments()),
                        internalNode.getFreeEnergyLB(bound.sequence),
                        internalNode.getFreeEnergyUB(bound.sequence)
                ));


                /*
                If the conformation has negative free energy and the bounds on that conf are loose enough that
                the subtree epsilon is worse than the targetEpsilon,
                then dive in a DFS manner to quickly tighten the node bound
                 */

                if (internalNode.getFreeEnergyUB(bound.sequence) < 0 &&
                        internalNode.getFreeEnergyLB(bound.sequence) - startingELB > drillDownDifference) {
                    //TODO: Put this back in possibly
                    /*
                    loopTasks.submit(() -> {
                        boundLowestBoundConfUnderNode(bound, internalNode, newNodes);
                        return null;
                    }, (ignored) -> {
                    });

                     */
                    System.out.println("Found drilldown opportunity");
                    processPartialConfNode(bound, newNodes, internalNode);
                } else {
                    processPartialConfNode(bound, newNodes, internalNode);
                }
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

        double endingELB = bound.calcEBound(e -> e.getFreeEnergyLB(bound.sequence));
        double endingEUB = bound.calcEBound(e -> e.getFreeEnergyUB(bound.sequence));
        double endingEps = bound.calcEpsilon();

        System.out.println(String.format("Ending error: %.10f, spread of [%.3f, %.3f]",
                endingEps,
                endingELB,
                endingEUB));

        if (endingEps > startingEps) {
            pilotFish.travelTree(bound.sequence);
            throw new RuntimeException("ERROR! bounds got looser");
        }
    }

    /**
     * Get two lists of nodes that we would like to process
     *
     * @param seqBound      The single-sequence pfunc bound
     * @param internalNodes The list of internal nodes to return
     * @param leafNodes     The list of leaf nodes to return
     * @param zSums         Tracking variables for the above lists containing the partition function error
     *                      Associated with each list
     *                      <p>
     *                      TODO: implement me
     */
    public void populateQueues(SingleSequenceSHARKStarBound_refactor seqBound, List<SHARKStarNode> internalNodes, List<SHARKStarNode> leafNodes, BigDecimal[] zSums) {
        List<SHARKStarNode> leftoverLeaves = new ArrayList<>();
        //PriorityQueue<SHARKStarNode> queue = seqBound.fringeNodes;
        // temporary fix
        PriorityQueue<SHARKStarNode> queue = null;

        //int maxNodes = 1000;
        int maxNodes = 10;
        if (leafTimeAverage > 0)
            maxNodes = Math.max(maxNodes, (int) Math.floor(0.1 * leafTimeAverage / internalTimeAverage));
        while (!queue.isEmpty() && (seqBound.internalQueue.size() < maxNodes
                || (!seqBound.leafQueue.isEmpty() && queue.peek().getScore(seqBound.sequence) >
                seqBound.leafQueue.peek().getScore(seqBound.sequence))
                || (!seqBound.internalQueue.isEmpty() && queue.peek().getScore(seqBound.sequence) >
                seqBound.internalQueue.peek().getScore(seqBound.sequence))
                || seqBound.leafQueue.size() < seqBound.maxMinimizations)) {
            SHARKStarNode curNode = queue.poll();

            // TODO: are these next two lines necessary at all?
            ConfIndex index = new ConfIndex(fullRCs.getNumPos());
            curNode.index(index);

            // Apply corrections to nodes
            /*
            boolean nodeCorrected = applyCorrectionsOrNOOP(curNode, seqBound);

            // If we corrected a node, then move to the next node
            if (nodeCorrected){
                leftoverLeaves.add(curNode);
                continue;
            }
             */

            if (curNode.getLevel() < fullRCs.getNumPos()) {
                seqBound.internalQueue.add(curNode);
                //} else if (shouldMinimize(curNode, seqBound.sequence) && !correctedNode(leftoverLeaves, curNode, node, seqBound.sequence)) {
                //TODO: Make sure the correctedNode method is doing what I think it is
                //} else if (shouldMinimize(curNode, seqBound.sequence) && !nodeCorrected){
            } else {
                seqBound.leafQueue.add(curNode);
            }

        }

        zSums[0] = fillListFromQueue(internalNodes, seqBound.internalQueue, maxNodes, seqBound);
        zSums[1] = fillListFromQueue(leafNodes, seqBound.leafQueue, seqBound.maxMinimizations, seqBound);
        if (!seqBound.internalQueue.isEmpty() &&
                MathTools.isLessThan(zSums[0], bc.calc(seqBound.internalQueue.peek().getScore(seqBound.sequence))))
            System.out.println("Should have used a node from the internal queue. How??");
        queue.addAll(leftoverLeaves);
    }

    private BigDecimal fillListFromQueue(List<SHARKStarNode> list, Queue<SHARKStarNode> queue, int max, SingleSequenceSHARKStarBound_refactor bound) {
        BigDecimal sum = BigDecimal.ZERO;
        List<SHARKStarNode> leftovers = new ArrayList<>();
        while (!queue.isEmpty() && list.size() < max) {
            SHARKStarNode curNode = queue.poll();
            /* change this to use the new applyCorrections method
            if (correctedNode(leftovers, curNode, curNode.getConfSearchNode(), seq)) {
                BigDecimal diff = curNode.getUpperBound(seq).subtract(curNode.getLowerBound(seq));
                if(MathTools.isGreaterThan(diff, sum)) {
                    leftovers.remove(curNode);
                    queue.add(curNode);
                }
                continue;
            }
             */
            /*
            if (applyCorrectionsOrNOOP(curNode, bound)){
                BigDecimal diff = curNode.getUpperBound(bound.sequence).subtract(curNode.getLowerBound(bound.sequence));
                if(MathTools.isGreaterThan(diff, sum)) {
                    leftovers.remove(curNode);
                    queue.add(curNode);
                }
                continue;
            }
             */
            BigDecimal diff = bc.calc(curNode.getFreeEnergyLB(bound.sequence)).subtract(bc.calc(curNode.getFreeEnergyUB(bound.sequence)));
            sum = sum.add(diff);
            list.add(curNode);
        }
        queue.addAll(leftovers);
        return sum;
    }

    /**
     * Do scoring and checks for a conformation with all residues assigned a rotamer
     *
     * @param seqBound     The single-sequence pfunc bound
     * @param newNodes     A list to track the new nodes created by this action
     * @param fullConfNode The full conformation node to process
     *                     <p>
     *                     TODO: implement me
     */
    public void processFullConfNode(SingleSequenceSHARKStarBound_refactor seqBound, List<SHARKStarNode> newNodes, SHARKStarNode fullConfNode) {
        //Sanity check, this should eventually never happen
        if (fullConfNode.getFreeEnergyLB(seqBound.sequence) > 10 &&
                (!seqBound.internalQueue.isEmpty() && seqBound.internalQueue.peek().getFreeEnergyLB(seqBound.sequence) < 0
                        || !seqBound.leafQueue.isEmpty() && seqBound.leafQueue.peek().getFreeEnergyLB(seqBound.sequence) < 0)) {
            System.err.println("not processing high-energy conformation");
            newNodes.add(fullConfNode);
            return;
        }

        //TODO: move this into its own method

        // Try to apply partial minimization correction to child
        boolean didCorrect = false;
        // Try LB
        double HOTCorrection = lowerBoundCorrector.getCorrection(fullConfNode.getAssignments());
        if (HOTCorrection > fullConfNode.getHOTCorrectionLB()){
            fullConfNode.setHOTCorrectionLB(HOTCorrection);
            didCorrect=true;
            // if we applied a correction, done with the node
        }
        if (didCorrect == true){
            newNodes.add(fullConfNode);
            return;
        }

        // Minimize the node
        minimizeNodeForSeq(fullConfNode, seqBound.sequence);
        //progress.reportLeafNode(fullConfNode.getMinE(), seqBound.fringeNodes.size(), seqBound.calcEpsilon());

        // Add back to leaf Queue
        //seqBound.leafQueue.add(fullConfNode);
        newNodes.add(fullConfNode);
    }

    /**
     * Do scoring and checks for a conformation with some residues unassigned
     *
     * @param seqBound        The single-sequence pfunc bound
     * @param newNodes        A list to track the new nodes created by this action
     * @param partialConfNode The partial conformation node to process
     */
    public void processPartialConfNode(SingleSequenceSHARKStarBound_refactor seqBound, List<SHARKStarNode> newNodes, SHARKStarNode partialConfNode) {
        //TODO: move this into its own method
        //Try just correcting it
        boolean didCorrect = false;
        // Try correcting LB
        double HOTCorrectionLB = lowerBoundCorrector.getCorrection(partialConfNode.getAssignments());
        if(HOTCorrectionLB > partialConfNode.getHOTCorrectionLB()){
            partialConfNode.setHOTCorrectionLB(HOTCorrectionLB);
            didCorrect = true;
        }
        if (didCorrect == true){
            newNodes.add(partialConfNode);
            return;
        }

        // Get the next position
        int nextPos;
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            partialConfNode.index(context.index);
            nextPos = order.getNextPos(context.index, seqBound.seqRCs);
        }

        // score child nodes with tasks (possibly in parallel)
        List<SHARKStarNode> children = new ArrayList<>();
        for (int nextRC : seqBound.seqRCs.get(nextPos)) {
            // Actually score the child
            makeAndScoreChildForSeq(partialConfNode, nextPos, nextRC, seqBound.sequence, seqBound.seqRCs, children);
        }
        loopTasks.waitForFinish();
        //children.stream().map(c -> c.toSeqString(seqBound.sequence)).forEach(System.out::println);
        partialConfNode.addChildren(children);
        // reporting
        //TODO: maybe separate scoring tasks from minimization tasks?
        loopTasks.waitForFinish();
        //TODO: fix this reporting process, it's almost certainly slowing things down
        /*
        for (SHARKStarNode child : children){
            progress.reportInternalNode(child.getLevel(), child.getPartialConfLB(),
                    child.getFreeEnergyLB(seqBound.sequence), seqBound.fringeNodes.size(), children.size(),
                    seqBound.calcEpsilon());
        }
         */
        newNodes.addAll(children);
    }

    /**
     * Expand partial conformation nodes until either the numNodesCeiling is met or the next partial Conf error is below the errorFloor
     * @param bound
     * @param numNodesCeiling   A soft ceiling on the number of nodes we are allowed to create this round
     * @param errorFloor
     */
    private PartialResult processPartialConfNodes(SingleSequenceSHARKStarBound_refactor bound, int numNodesCeiling, double errorFloor) {
        PartialResult result = new PartialResult();
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();

            Stopwatch partialWatch = new Stopwatch().start();

            int numNodesCreated = 0;
            while (numNodesCreated < numNodesCeiling && !bound.internalQueue.isEmpty()) {     // Stop if we create enough nodes, or if there are no more nodes
                SHARKStarNode parent = bound.internalQueue.take();

                // Stop if the node error falls below our cutoff
                if(parent.getScore(bound.sequence) < errorFloor){
                    bound.internalQueue.put(parent);
                    return result;
                }

                //Record starting parent energy bounds
                double parentLB = parent.getFreeEnergyLB(bound.sequence);
                double parentUB = parent.getFreeEnergyUB(bound.sequence);

                // Try to correct node
                boolean correctedParent = correctNodeOrNOOP(parent, parent.getHOTCorrectionLB(), bound, context);
                if(correctedParent){
                    result.lowerBoundDeltas.add(bc.calc_EDiff(parentLB, parent.getFreeEnergyLB(bound.sequence)));
                    result.upperBoundDeltas.add(bc.calc_EDiff(parentUB, parent.getFreeEnergyUB(bound.sequence)));
                    bound.internalQueue.put(parent);
                    numNodesCreated++; // technically we didn't create it, but we updated
                    continue;
                }

                parent.index(context.index);
                int nextPos = order.getNextPos(context.index, bound.seqRCs);

                // Create list for finding out the delta pfunc upper and lower
                ArrayList<Double> childLBs = new ArrayList<>();
                ArrayList<Double> childUBs = new ArrayList<>();

                // score child nodes
                int[] nextRCs = bound.seqRCs.get(nextPos);
                for (int nextRC : nextRCs) {
                    // Make the child node
                    parent.index(context.index);    // this is because the index will change in the loop
                    SHARKStarNode child = parent.assign(nextPos, nextRC);
                    // Score the child node
                    child.setPartialConfLB(context.partialConfLBScorer.calcDifferential(context.index, bound.seqRCs, nextPos, nextRC));
                    child.setPartialConfUB(context.partialConfUBScorer.calcDifferential(context.index, bound.seqRCs, nextPos, nextRC));
                    child.index(context.index); // We don't implement calcDifferential for the SHARKStarNodeScorer yet, so probably faster to do this
                    child.setUnassignedConfLB(context.unassignedConfLBScorer.calc(context.index, bound.seqRCs), bound.sequence);
                    child.setUnassignedConfUB(context.unassignedConfUBScorer.calc(context.index, bound.seqRCs), bound.sequence);

                    // Try to apply corrections to child
                    correctNodeOrNOOP(child, parent.getHOTCorrectionLB(), bound, context);

                    // Now we have good bounds on child energy
                    double childLB = child.getFreeEnergyLB(bound.sequence);
                    double childUB = child.getFreeEnergyUB(bound.sequence);

                    // Compute the node partition function error
                    child.setScore(
                            bc.calc_lnZDiff(childLB, childUB),
                            bound.sequence);
                    // Record child bounds for use in this method
                    childLBs.add(childLB);
                    childUBs.add(childUB);

                    // Add child to parent's children ( we can get rid of this eventually)
                    parent.addChild(child);

                    // Add child back to queue
                    if(child.getLevel() < bound.seqRCs.getNumPos())
                        bound.internalQueue.put(child);
                    else{
                        bound.leafQueue.put(child);
                    }

                    numNodesCreated++;
                }

                //Now we are done making children for the parent node, let's update the result
                //Compute the delta E bounds
                result.lowerBoundDeltas.add(bc.calc_EDiff(parentLB, bc.logSumExp(childLBs)));
                result.upperBoundDeltas.add(bc.calc_EDiff(bc.logSumExp(childUBs), parentUB));
            }
            partialWatch.stop();
            result.time = partialWatch.getTimeS();
            result.numNodesCreated = numNodesCreated;
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return result;
    }

    public MinResult minimizeNode(SingleSequenceSHARKStarBound_refactor bound){
        MinResult result = new MinResult();
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            //Timing
            SHARKStarNode node = bound.leafQueue.take();
            Stopwatch minimizationTimer = new Stopwatch();
            Stopwatch correctionTimer = new Stopwatch();

            ScoreContext context = checkout.get();
            node.index(context.index);

            double prevLB = node.getFreeEnergyLB(bound.sequence);
            double prevUB = node.getFreeEnergyUB(bound.sequence);

            ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.getAssignments(), node.getFreeEnergyLB(bound.sequence));
            minimizationTimer.start();
            ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);
            node.setMinE(analysis.epmol.energy);
            node.setIsMinimized(true);
            node.setScore(Double.NEGATIVE_INFINITY, bound.sequence);

            System.out.println(String.format("Minimized %s --> %.3f in [%.3f, %.3f]",
                    node.confToString(),
                    node.getMinE(),
                    prevLB,
                    prevUB
            ));

            // Compute deltas
            result.lowerBoundDelta = bc.calc_EDiff(prevLB, analysis.epmol.energy);
            result.upperBoundDelta = bc.calc_EDiff(analysis.epmol.energy, prevUB);

            if (this.saveEPMOLsForMinimization) {
                synchronized (this) {
                    if (analysis.epmol != null) {
                        epmolsMap.put(node.getAssignments(), analysis.epmol);
                    }
                }
            }

            minimizationTimer.stop();
            bound.leafQueue.put(node);
            if (doExtraTupleCorrections) {
                correctionTimer.start();
                energyMatrixCorrector.computeEnergyCorrection(analysis, conf, 1.0,
                        context.batcher);
                correctionTimer.stop();
            }
            result.minimizationTime = minimizationTimer.getTimeS();
            result.correctionTime = correctionTimer.getTimeS();


        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return result;
    }

    /**
     * Try to apply corrections to node.
     *
     * SIDE EFFECT: Changes context.index
     *
     * @param node
     * @param lowerCorrectionFloor
     * @param bound
     * @param context
     * @return
     */
    private boolean correctNodeOrNOOP(SHARKStarNode node, double lowerCorrectionFloor, SingleSequenceSHARKStarBound_refactor bound, ScoreContext context){
        boolean corrected = false;
        node.index(context.index);
        // Try to apply partial minimization correction to child (Lowerbound)
        // If the parent correction is larger, just use the parent correction, since it must be valid for lowerbound
        double HOTCorrectionLB = Math.max(lowerBoundCorrector.getCorrection(context.index), lowerCorrectionFloor);
        if (HOTCorrectionLB > node.getHOTCorrectionLB()){
            node.setHOTCorrectionLB(HOTCorrectionLB);
            corrected = true;
        }

        // Try to apply correction to child (upperbound)
        if(doUpperBoundCorrections && !saveEPMOLsForMinimization && upperBoundCorrector.getEntries(new RCTuple(context.index.makeConf())).size()>0) {
            RCTuple mappedTup = upperBoundCorrector.getEntries(new RCTuple(context.index.makeConf())).get(0).mapTup;
            mappedTup.pasteToIndex(context.index);
            double correctedPartialUB = context.upperCorrectionGScorer.calc(context.index, bound.seqRCs);
            double correctedUnassignedUB = context.upperCorrectionHScorer.calc(context.index, bound.seqRCs);
            double HOTCorrectionUB = correctedPartialUB + correctedUnassignedUB - node.getFreeEnergyUB(bound.sequence);
            // Filter out positive upper bound corrections
            //TODO: make this sequence dependent
            if (HOTCorrectionUB < 0) {
                node.setHOTCorrectionUB(HOTCorrectionUB);
                corrected = true;
            }
        }
        node.index(context.index);
        return corrected;
    }

    /**
     * Drill down in a standard A* manner and bound the lowest-energy full conformation in the subtree
     * defined by internalNode
     *
     * @param seqBound     The single-sequence pfunc bound
     * @param internalNode The node defining the subtree
     * @param newNodes     A list to which we add the new nodes generated by this process
     */
    public void boundLowestBoundConfUnderNode(SingleSequenceSHARKStarBound_refactor seqBound, SHARKStarNode internalNode, List<SHARKStarNode> newNodes) {
        System.out.println("Bounding " + internalNode.toSeqString(seqBound.sequence));

        // Set up custom queue that mimics usual A* behavior
        Comparator<SHARKStarNode> confBoundComparator = Comparator.comparingDouble(o -> o.getFreeEnergyLB(seqBound.sequence));
        PriorityQueue<SHARKStarNode> drillQueue = new PriorityQueue<>(confBoundComparator);
        drillQueue.add(internalNode);

        List<SHARKStarNode> generatedNodes = new ArrayList<>();
        int numNodes = 0;
        Stopwatch leafLoop = new Stopwatch().start();
        Stopwatch overallLoop = new Stopwatch().start();
        while (!drillQueue.isEmpty()) {
            numNodes++;
            SHARKStarNode curNode = drillQueue.poll();
            ConfIndex index = new ConfIndex(seqBound.seqRCs.getNumPos());
            curNode.index(index);

            if (curNode.getLevel() < seqBound.seqRCs.getNumPos()) {
                SHARKStarNode nextNode = drillDown(seqBound, generatedNodes, curNode);
                // Sometimes there are no good leaf nodes. Weird.
                if (nextNode != null) {
                    generatedNodes.remove(nextNode);
                    drillQueue.add(nextNode);
                }
            } else {
                generatedNodes.add(curNode);
            }

            //debugHeap(drillQueue, true);
            if (leafLoop.getTimeS() > 1) {
                leafLoop.stop();
                leafLoop.reset();
                leafLoop.start();
                System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",
                        numNodes,
                        overallLoop.getTime(2),
                        seqBound.calcZBound(e -> e.getFreeEnergyUB(seqBound.sequence)),
                        seqBound.calcZBound(e -> e.getFreeEnergyLB(seqBound.sequence))
                ));
            }
        }
        newNodes.addAll(generatedNodes);
    }

    /**
     * A* search using energy lower bound for a leaf node
     *
     * @param seqBound       The single-sequence pfunc bound
     * @param generatedNodes A list to which we add new nodes found during this process
     * @param node           The node to start the drillDown process from
     * @return the lowest-energy leaf node
     * <p>
     * TODO: Implement me
     */
    public SHARKStarNode drillDown(SingleSequenceSHARKStarBound_refactor seqBound, List<SHARKStarNode> generatedNodes, SHARKStarNode node) {
        throw new NotImplementedException();
    }

    /**
     * Do cleanup after one loop iteration
     *
     * @param seqBound The single-sequence pfunc bound
     * @param newNodes A list of new nodes created during the loop iteration
     * @param timer    The timer telling us how long the loop took
     * @param numNodes The number of nodes we processed during the loop
     *                 <p>
     *                 TODO: implement me
     */
    public void loopCleanup(SingleSequenceSHARKStarBound_refactor seqBound, List<SHARKStarNode> newNodes, Stopwatch timer, int numNodes) {
        //PriorityQueue<SHARKStarNode> queue = seqBound.fringeNodes;
        // Temporary fix
        PriorityQueue<SHARKStarNode> queue = null;
        for (SHARKStarNode node : newNodes) {
            if (node != null) {
                if (node.isMinimized())
                    seqBound.addFinishedNode(node);
                else
                    queue.add(node);
            }
        }
        timer.stop();
        double loopTime = timer.getTimeS();
        System.out.println("Processed " + numNodes + " this loop, spawning " + newNodes.size() + " in " + loopTime);
        timer.reset();
        timer.start();
        //energyMatrixCorrector.processPreminimization(seqBound, minimizingEcalc);
        //profilePrint("Preminimization time : " + timer.getTime(2));
        //printTree(seqBound.sequence, rootNode);
        timer.stop();
        //cleanupTime = timer.getTimeS();
        //double scoreChange = rootNode.updateAndReportConfBoundChange(new ConfIndex(RCs.getNumPos()), RCs, correctiongscorer, correctionhscorer);

        System.out.println(String.format("Loop complete. Epsilon is now %.10f, Bounds are now [%.3f,%.3f]", seqBound.calcEpsilon(), seqBound.calcEBound(e -> e.getFreeEnergyLB(seqBound.sequence)),
                seqBound.calcEBound(e -> e.getFreeEnergyUB(seqBound.sequence))));
    }

    protected void debugPrint(String s) {
        if (debug)
            System.out.println(s);
    }

    public void printEnsembleAnalysis() {
        ensembleAnalyzer.printStats();
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

    public EnergyMatrix getMinimizingEmat() {
        return minimizingEmat;
    }

    public ConfEnergyCalculator getMinimizingEcalc() {
        return minimizingEcalc;
    }

    public Set<String> getCorrectedTuples() {
        return correctedTuples;
    }

    public SimpleScoreCorrector getCorrectionMatrix() {
        return lowerBoundCorrector;
    }

    public void setNumPartialMinimizations(int numPartialMinimizations) {
        this.numPartialMinimizations = numPartialMinimizations;
    }

    public void setComputedCorrections(boolean computedCorrections) {
        this.computedCorrections = computedCorrections;
    }

    protected static class ScoreContext {
        public ConfIndex index;
        public AStarScorer partialConfLBScorer;
        public AStarScorer partialConfUBScorer;
        public AStarScorer unassignedConfLBScorer;
        public AStarScorer unassignedConfUBScorer;
        public ConfEnergyCalculator ecalc;
        public BatchCorrectionMinimizer batcher;
        public AStarScorer upperCorrectionGScorer;
        public AStarScorer upperCorrectionHScorer;
    }

    private class ScoringResult {
        SHARKStarNode resultNode = null;
        double partialLB = Double.NaN;
        double partialUB = Double.NaN;
        double unassignLB = Double.NaN;
        double unassignUB = Double.NaN;
        double HOTCorrectionLB = 0.0;
        double HOTCorrectionUB = 0.0;
        double score = Double.NaN;
        double time = Double.NaN;
        String historyString = "Error!!";

        public boolean isValid() {
            return resultNode != null && !Double.isNaN(partialLB) && !Double.isNaN(partialUB) && !Double.isNaN(unassignLB) && !Double.isNaN(unassignUB) && !Double.isNaN(score);
        }
    }

    private class MinimizationResult {
        SHARKStarNode resultNode = null;
        double minimizedEnergy = Double.NaN;
        double minimizationTime = Double.NaN;
        double correctionTime = Double.NaN;
        EnergyCalculator.EnergiedParametricMolecule epmol = null;
        String historyString = "Error!!";

        public boolean isValid() {
            return resultNode != null && !Double.isNaN(minimizedEnergy);
        }
    }

    private class PartialResult{
        double time;
        int numNodesCreated;
        ArrayList<Double> lowerBoundDeltas = new ArrayList<>();
        ArrayList<Double> upperBoundDeltas = new ArrayList<>();
    }

    private class MinResult{
        double minimizationTime;
        double correctionTime;
        int numNodesCreated = 1;
        double lowerBoundDelta;
        double upperBoundDelta;
    }

    public static interface RigidEmatFactory{
        EnergyMatrix make(SimpleConfSpace confSpace);
    }

    private static enum Step {
        None,
        Score,
        Energy
    }
}

