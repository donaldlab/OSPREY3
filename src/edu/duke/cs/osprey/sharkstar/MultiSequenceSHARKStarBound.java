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

import java.io.BufferedWriter;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.printTree;

public class MultiSequenceSHARKStarBound implements PartitionFunction {

    private final EnergyMatrixCorrector energyMatrixCorrector;
    private Sequence precomputedSequence;
    protected double targetEpsilon = 1;
    public static final boolean debug = false;
    public static final boolean suppressPrecisionWarnings = true;
    public static final boolean diveForLeaves = false;
    public boolean profileOutput = false;
    private Status status = null;

    // the number of full conformations minimized
    private int numConfsEnergied = 0;
    private int numConfsEnergiedThisLoop = 0;

    // the number of full conformations scored OR energied
    private int numConfsScored = 0;
    // max num confs to minimize, -1 is unbounded
    private int maxNumConfs = -1;
    private boolean nonZeroLower;

    protected int numInternalNodesProcessed = 0;
    private boolean computedCorrections = false;

    private boolean printMinimizedConfs = false;
    private final MARKStarProgress progress;
    public String stateName = String.format("%4f", Math.random());
    private int numPartialMinimizations;
    public ArrayList<Integer> minList;
    protected double internalTimeAverage;
    protected double leafTimeAverage;
    protected static TaskExecutor loopTasks;


    // We keep track of the root node for computing our K* bounds
    public MultiSequenceSHARKStarNode rootNode;
    // Heap of nodes for recursive expansion
    private final ConfIndex confIndex;
    public StaticBiggestLowerboundDifferenceOrder order;
    public final AStarPruner pruner;
    // TODO: Implement new AStarPruner for MARK*?
    protected RCs fullRCs;
    protected Parallelism parallelism;
    public ObjectPool<ScoreContext> contexts;
    private final ScorerFactory gscorerFactory;
    private final ScorerFactory rigidgscorerFactory;
    private final ScorerFactory hscorerFactory;
    private final ScorerFactory nhscorerFactory;

    private final ConfAnalyzer confAnalyzer;
    EnergyMatrix minimizingEmat;
    EnergyMatrix rigidEmat;
    UpdatingEnergyMatrix correctionMatrix;
    ConfEnergyCalculator minimizingEcalc;
    final BatchCorrectionMinimizer theBatcher;

    private final SHARKStarEnsembleAnalyzer ensembleAnalyzer;
    private final Stopwatch stopwatch = new Stopwatch().start();
    // Variables for reporting pfunc reductions more accurately
    BigDecimal startUpperBound;
    BigDecimal startLowerBound;
    BigDecimal upperReduction_PartialMin = BigDecimal.ZERO; //Pfunc upper bound improvement from partial minimization corrections

    BigDecimal cumulativeZCorrection = BigDecimal.ZERO;//Pfunc upper bound improvement from partial minimization corrections
    BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
    private final Set<String> correctedTuples = Collections.synchronizedSet(new HashSet<>());
    private BigDecimal stabilityThreshold;
    private double internalTimeSum = 0;
    private final MultiSequenceState state;


    private MultiSequenceSHARKStarBound precomputedPfunc;
    public MultiSequenceSHARKStarNode precomputedRootNode;
    public final SimpleConfSpace confSpace;

    private BigDecimal precomputedUpperBound;
    private BigDecimal precomputedLowerBound;

    private final List<MultiSequenceSHARKStarNode> precomputedFringe = new ArrayList<>();

    public static final int[] debugConf = new int[]{};//6, 5, 15, -1, -1, 8, -1};
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

        rootNode.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, bc.calc(confLowerBound), bc.calc(confUpperBound), precomputedSequence, "(root initialization)");

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
        //computeFringeForSequence(newBound, this.rootNode);

        // Compute the fringe in parallel
        List<MultiSequenceSHARKStarNode> scoredFringe = computeFringeForSequenceParallel(newBound.sequence, newBound.seqRCs);
        if(scoredFringe.size()==0)
            scoredFringe.add(this.rootNode);
        debugPrint(String.format("[Normal fringe # nodes, Parallel fringe # nodes] = [%d, %d]",newBound.fringeNodes.size(), scoredFringe.size()));
        newBound.fringeNodes.addAll(scoredFringe);

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
        if (node.hasChildren(precomputedSequence)) {
            for (MultiSequenceSHARKStarNode child : node.getChildren(precomputedSequence)) {
                updatePrecomputedNode(child, permutation, size);
            }
        }else {
            precomputedFringe.add(node);
        }
    }

    private static class FringeResult{
        boolean isFringe;
        List<MultiSequenceSHARKStarNode> nodes = new ArrayList<>();
    }

    private List<MultiSequenceSHARKStarNode> computeFringeForSequenceParallel(Sequence seq, RCs seqRCs) {
        List<MultiSequenceSHARKStarNode> scoredSeqFringe = Collections.synchronizedList(new ArrayList<>());
        ConcurrentLinkedQueue<MultiSequenceSHARKStarNode> unscoredSeqFringe = new ConcurrentLinkedQueue<>(precomputedFringe);

        // Sometimes the precomputedFringe will be empty. It really shouldn't be, but for now just add the root Node if it is
        if(unscoredSeqFringe.isEmpty())
            unscoredSeqFringe.add(this.rootNode);

        while(!unscoredSeqFringe.isEmpty() || loopTasks.isExpecting()){
            MultiSequenceSHARKStarNode node = unscoredSeqFringe.poll();
            if(node == null)
                continue;

            loopTasks.submitExpecting(
                ()->{
                    FringeResult result = new FringeResult();
                    try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                        ScoreContext context = checkout.get();

                        // Fix issue where nextDesignPosition can be null
                        MultiSequenceSHARKStarNode.Node confNode = node.getConfSearchNode();
                        confNode.index(context.index);
                        if (node.nextDesignPosition == null && node.level < confSpace.positions.size()) {
                            node.nextDesignPosition = confSpace.positions.get(order.getNextPos(context.index, seqRCs));
                        }

                        // Get the children
                        List<MultiSequenceSHARKStarNode> children = node.getChildren(seq);
                        // If we are at a leaf, score the node
                        if(children != null && !children.isEmpty()) {
                            // if there are children, just add them to queue, since we only want the fringe
                            result.isFringe = false;
                            result.nodes.addAll(children);
                            //unscoredSeqFringe.addAll(children);
                        }else{
                            double confCorrection = correctionMatrix.confE(confNode.assignments);
                            double gscore = context.partialConfLowerBoundScorer.calc(context.index, seqRCs);
                            double hscore = context.lowerBoundScorer.calc(context.index, seqRCs);
                            double confLowerBound = confNode.getPartialConfLowerBound() + context.lowerBoundScorer.calc(context.index, seqRCs);
                            double confUpperBound = confNode.getPartialConfUpperBound() + context.upperBoundScorer.calc(context.index, seqRCs);
                            String historyString = "";
                            if(debug){
                                historyString = String.format("%s: previous lower bound %f, g score %f, hscore %f, f score %f corrected score %f, from %s",
                                        confNode.confToString(), node.getConfLowerBound(seq), gscore, hscore, gscore + hscore, confCorrection, getStackTrace());
                            }
                            node.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound,bc.calc(confLowerBound), bc.calc(confUpperBound),  seq, historyString);

                            if (node.getChildren(null).isEmpty())
                                correctionMatrix.setHigherOrder(node.toTuple(), confNode.getPartialConfLowerBound()
                                        - minimizingEmat.confE(confNode.assignments));

                            result.isFringe = true;
                            result.nodes.add(node);
                            //bound.fringeNodes.add(node);
                            //scoredSeqFringe.add(node);
                        }
                    }
                    return result;
                },
                (FringeResult result)->{
                    if (result.isFringe)
                        scoredSeqFringe.addAll(result.nodes);
                    else
                        unscoredSeqFringe.addAll(result.nodes);
                }
            );
        }
        loopTasks.waitForFinish();
        return scoredSeqFringe;

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

    private class CorrectionResult{
       boolean didCorrect = false;
       BigDecimal deltaUB;
       double correctionSize;
       MultiSequenceSHARKStarNode correctedNode;
    }

    /**
     * Try to apply corrections to node.
     */
    private CorrectionResult correctNodeOrFalse(MultiSequenceSHARKStarNode node, SingleSequenceSHARKStarBound bound) {
        CorrectionResult result = new CorrectionResult();
        double confCorrection = correctionMatrix.confE(node.getConfSearchNode().assignments);
        double oldg = node.getConfSearchNode().getPartialConfLowerBound();
        double correctionDiff = confCorrection - oldg;

        if ( correctionDiff > 1e-5) {
            result.didCorrect = true;

            BigDecimal oldZUpperBound = node.getUpperBound(bound.sequence);
            double oldConfLowerBound = node.getConfLowerBound(bound.sequence);

            // update the node gscore
            node.getConfSearchNode().setPartialConfLowerAndUpper(confCorrection, node.getConfSearchNode().getPartialConfUpperBound());
            recordCorrection(oldg, correctionDiff);
            String historyString = String.format("%s: correction from %f to %f, from ",
                    node.getConfSearchNode().confToString(), oldg, confCorrection);

            // update the node total scores
            node.setBoundsFromConfLowerAndUpperWithHistory(oldConfLowerBound + correctionDiff, node.getConfUpperBound(bound.sequence), bc.calc(oldConfLowerBound+correctionDiff), bc.calc(node.getConfUpperBound(bound.sequence)), bound.sequence, historyString);
            node.markUpdated();
            debugPrint("Correcting " + node.toSeqString(bound.sequence) +" correction ="+(correctionDiff) );

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

            result.deltaUB = diffUB;
            result.correctedNode = node;
        }
        return result;
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
        System.out.println("Tightening bound for "+sequenceBound.sequence + " " + this.stateName);

        // Test to make sure that we have a valid pfunc
        sequenceBound.updateBound();
        double lastEps = sequenceBound.getSequenceEpsilon();
        if(lastEps == 0)
            System.err.println("Epsilon is ZERO??! And we are still tightening the bound!?");

        int previousConfCount = workDone();

        // Run to populate the queues a little bit
        if (!sequenceBound.nonZeroLower()) {
            runUntilNonZero(sequenceBound);
            sequenceBound.updateBound();
        }

        // Initialize state

        /*
        sequenceBound.state.upperBound = sequenceBound.fringeNodes.stream().map(n -> n.getUpperBound(sequenceBound.sequence)).
                reduce(BigDecimal.ZERO, (a,b) -> a.add(b, PartitionFunction.decimalPrecision));
        sequenceBound.state.lowerBound = sequenceBound.fringeNodes.stream().map(n -> n.getLowerBound(sequenceBound.sequence)).
                reduce(BigDecimal.ZERO, (a,b) -> a.add(b, PartitionFunction.decimalPrecision));

         */

        sequenceBound.state.upperBound = sequenceBound.getUpperBound();
        sequenceBound.state.lowerBound = sequenceBound.getLowerBound();
        double curEps = sequenceBound.state.calcDelta();

        Step step = Step.None;

        /*
        Continue looping if 1) there are nodes in the fringe queue or
        2) we are running tasks that will result in nodes being added to the queue
         */
        this.numConfsEnergiedThisLoop = 0;

        int numIters = 0;
        while( (!sequenceBound.fringeNodes.isEmpty() || loopTasks.isExpecting())){
            numIters++;
            if(numIters> 10000){
                break;
            }

            synchronized(sequenceBound) {
                //TODO: apply partial minimizations
                sequenceBound.updateBound();

                // Early termination
                double newEps = sequenceBound.state.calcDelta();
                if(newEps > curEps){
                    //throw new RuntimeException("ERROR: Epsilon is increasing");
                }
                curEps = newEps;

                debugPrint(String.format("Epsilon: %.9f, Bounds:[%1.3e, %1.3e]",
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
                        debugPrint("Exiting loop because of work");

                    if(!isStable(stabilityThreshold, sequenceBound)) {
                        debugPrint("Exiting loop due to stablity, thresh: " + stabilityThreshold);
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
            /*
            synchronized(sequenceBound){
                double leafTimeAvg = sequenceBound.state.totalTimeEnergy / sequenceBound.state.numEnergiedConfs;
                double internalTimeAvg = sequenceBound.state.totalTimeExpansion / sequenceBound.state.numExpansions;
                if (leafTimeAverage > 0)
                    maxInternalNodes = Math.min(maxInternalNodes, (int) Math.floor(0.1 * leafTimeAvg / internalTimeAvg));
                maxInternalNodes = Math.min(maxInternalNodes,
                        (int) (sequenceBound.fringeNodes.size() / 2*loopTasks.getParallelism()));
                maxInternalNodes = Math.max(maxInternalNodes,1);
            }
             */

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
                    MultiSequenceSHARKStarNode node = sequenceBound.fringeNodes.poll();

                    /*
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
                     */
                    // now, decide which step to take
                    /*
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
                     */
                    if(node.getConfSearchNode().getLevel() >= fullRCs.getNumPos()){
                        if(node.getErrorBound(sequenceBound.sequence).compareTo(BigDecimal.ONE) < 0){
                            System.err.println("Out of loose nodes!?");
                            sequenceBound.fringeNodes.add(node);
                            continue;
                        }else {
                            loosestLeaf = node;
                            step = Step.Energy;
                        }
                    }else{
                        internalNodes.add(node);
                        step = Step.ExpandInBatches;
                    }
                }
            }


            // Take steps
            switch(step) {
                case Energy:{
                    MultiSequenceSHARKStarNode toMinimize = loosestLeaf;
                    debugPrint(String.format("Minimizing node with lower bound %f and pfunc error %1.3e",
                            loosestLeaf.getConfLowerBound(sequenceBound.sequence),
                            loosestLeaf.getErrorBound(sequenceBound.sequence)
                            ));

                    loopTasks.submit( () -> {
                        CorrectionResult correctionResult = correctNodeOrFalse(toMinimize, sequenceBound);
                        if(correctionResult.didCorrect)
                            return correctionResult;
                        else
                            return minimizeNode(toMinimize, sequenceBound.sequence, sequenceBound.seqRCs);
                            },
                            (result) -> {
                        if (result.getClass() == CorrectionResult.class)
                            onCorrection((CorrectionResult) result, sequenceBound);
                        else
                            onMinimization((MinimizationResult) result, sequenceBound);
                            }
                    );
                    break;
                }
                case ExpandInBatches: {
                    List<MultiSequenceSHARKStarNode> toExpand = internalNodes;
                    numNodes = toExpand.size();
                    debugPrint("Processing " + numNodes + " internal nodes...");

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
                                    if (diveForLeaves && !MathTools.isGreaterThan(internalNode.getLowerBound(sequenceBound.sequence), BigDecimal.ONE) &&
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

                                BigDecimal lbAccuracyCutoff = startLB.multiply(BigDecimal.valueOf(-1e-10));
                                BigDecimal ubAccuracyCutoff = startUB.multiply(BigDecimal.valueOf(1e-10));

                                BigDecimal lbProportion = BigDecimal.ZERO;
                                BigDecimal ubProportion = BigDecimal.ZERO;
                                if(startLB.compareTo(BigDecimal.ZERO) > 0 && startUB.compareTo(BigDecimal.ZERO) > 0) {
                                    lbProportion = result.deltaLB.divide(startLB, PartitionFunction.decimalPrecision);
                                    ubProportion = result.deltaUB.divide(startUB, PartitionFunction.decimalPrecision);
                                }


                                if (result.deltaLB.compareTo(BigDecimal.ZERO) < 0) {
                                    //if the lower bound is decreasing at all
                                    //log it or something

                                    // If the wrong move magnitude is more than one
                                    if (result.deltaLB.compareTo(BigDecimal.valueOf(-1.0)) < 0) {

                                        // If the wrong move magnitude is more than a factor of 1e-10 of the original bound (and more than one), throw an exception
                                        if (result.deltaLB.compareTo(lbAccuracyCutoff) < 0) {
                                            System.out.println(String.format("ERROR: Expansion of %s resulted in increased UB\n\tdeltaLB: %1.9e, Starting LB: %1.9e, cutoff: %1.9e",
                                                    toExpand.get(0).toString(),
                                                    result.deltaLB,
                                                    startLB,
                                                    lbAccuracyCutoff
                                            ));
                                            throw new RuntimeException("Lower bound is decreasing");
                                        // If the wrong move magnitude is more than one but less than the cutoff, just warn
                                        } else {
                                            if(!suppressPrecisionWarnings) {
                                                System.err.println(String.format("WARNING: Expansion of %s resulted in LB decrease of %1.9e, a factor of %1.3e of the starting bound. This is likely just a numerical precision issue",
                                                        toExpand.get(0).toString(),
                                                        result.deltaLB,
                                                        lbProportion
                                                ));
                                            }
                                        }

                                    }
                                }
                                if (result.deltaUB.compareTo(BigDecimal.ZERO) > 0) {
                                    // if upper bound is increasing at all
                                    // log or something here

                                    // if the upper bound is increasing by magnitude greater than one
                                    if (result.deltaUB.compareTo(BigDecimal.ONE) > 0) {

                                        // If the wrong move magnitude is more than a factor of 1e-10 of the original bound (and more than one), throw an exception
                                        if (result.deltaUB.compareTo(ubAccuracyCutoff) > 0) {
                                            System.out.println(String.format("ERROR: Expansion of %s resulted in increased UB \n\tdeltaUB: %1.9e, Starting UB: %1.9e, cutoff: %1.9e",
                                                    toExpand.get(0).toString(),
                                                    result.deltaUB,
                                                    startUB,
                                                    ubAccuracyCutoff
                                            ));
                                            throw new RuntimeException("Upper bound is increasing");
                                        // If the wrong move magnitude is more than one but less than the cutoff, just warn
                                        }else{
                                            if(!suppressPrecisionWarnings) {
                                                System.err.println(String.format("WARNING: Expansion of %s resulted in UB increase of %1.9e, a factor of %1.3e of the starting bound. This is likely just a numerical precision issue",
                                                        toExpand.get(0).toString(),
                                                        result.deltaUB,
                                                        ubProportion
                                                ));
                                            }
                                        }
                                    }
                                }

                                internalTime.stop();
                                result.timeS = internalTime.getTimeS();
                                result.numExpanded = toExpand.size();
                                return result;
                            },
                            (result) -> {
                                double delta = 2.0;
                                long fringeSize = 0;
                                synchronized(sequenceBound.state) {
                                    // Update partition function values
                                    sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
                                    sequenceBound.state.lowerBound = sequenceBound.state.lowerBound.add(result.deltaLB, PartitionFunction.decimalPrecision);
                                    sequenceBound.state.numExpansions += result.numExpanded;
                                    sequenceBound.state.totalTimeExpansion += result.timeS;
                                    sequenceBound.state.numRoundsExpand++;
                                    delta = sequenceBound.state.calcDelta();
                                    fringeSize = sequenceBound.fringeNodes.size();
                                }

                                synchronized (this) {
                                    this.state.numExpansions += result.numExpanded;
                                    this.state.totalTimeExpansion += result.timeS;
                                    this.state.numRoundsExpand++;
                                    debugPrint("Got " + result.newNodes.size() +" internal nodes.");
                                    internalTimeSum = internalTime.getTimeS();
                                    internalTimeAverage = internalTimeSum / Math.max(1, toExpand.size());
                                    debugPrint("Internal node time :" + internalTimeSum + ", average " + internalTimeAverage);
                                    numInternalNodesProcessed += internalNodes.size();
                                    for (int i = 0; i< result.newNodes.size(); i++){
                                        MultiSequenceSHARKStarNode node = result.newNodes.get(i);
                                        this.progress.reportInternalNode(node.getConfSearchNode().getLevel(),
                                                node.getConfSearchNode().getGScore(),
                                                node.getConfLowerBound(sequenceBound.sequence) - node.getConfSearchNode().getGScore(),
                                                fringeSize,
                                                result.newNodes.size(),
                                                delta);
                                    }
                                    synchronized (sequenceBound.fringeNodes) {
                                        sequenceBound.fringeNodes.addAll(result.newNodes);

                                    }
                                }
                            }
                    );

                    break;
                }
                case Partial: {
                    debugPrint("Computing partial mins");
                    BatchCorrectionMinimizer.Batch batch = theBatcher.acquireBatch();
                    if(batch == null)
                        break;
                    loopTasks.submit(
                            () -> {
                                PartialMinimizationResult result = new PartialMinimizationResult();
                                Stopwatch partialMinTime = new Stopwatch().start();

                                // calculate all the fragment energies
                                Map<BatchCorrectionMinimizer.PartialMinimizationTuple, EnergyCalculator.EnergiedParametricMolecule> confs = new HashMap<>();
                                for (int i =0; i< batch.fragments.size(); i++) {
                                    BatchCorrectionMinimizer.PartialMinimizationTuple frag = batch.fragments.get(i);
                                    double energy;

                                    // are there any RCs are from two different backbone states that can't connect?
                                    if (theBatcher.isParametricallyIncompatible(frag.tup)) {

                                        // yup, give this frag an infinite energy so we never choose it
                                        energy = Double.POSITIVE_INFINITY;

                                    } else {

                                        // nope, calculate the usual fragment energy
                                        confs.put(frag, theBatcher.confEcalc.calcEnergy(frag.tup));
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
                                for(BatchCorrectionMinimizer.PartialMinimizationTuple tuple : result.confs.keySet()) {
                                    double lowerbound = theBatcher.minimizingEnergyMatrix.getInternalEnergy(tuple.tup);
                                    double uncorrectedParentConfLB = theBatcher.minimizingEnergyMatrix.getInternalEnergy(tuple.parentConf);
                                    double tupleEnergy = result.confs.get(tuple).energy;
                                    double correction = tupleEnergy - lowerbound;
                                    // Check to make sure that the tuple correction is not a minimizer error
                                    if (uncorrectedParentConfLB + correction - tuple.E > 1e-10) {
                                        System.err.println(String.format("WARNING: minimizer bug:\n\t%s has a correction of %f, but this would cause \n\tparent conf (%s) LB: %f (uncorrected %f) -> %f > parent conf E: %f",
                                                tuple.tup.toString(),
                                                correction,
                                                tuple.parentConf.toString(),
                                                tuple.parentConfLB,
                                                uncorrectedParentConfLB,
                                                correction+uncorrectedParentConfLB,
                                                tuple.E
                                                ));
                                    // Check to make sure that the tuple correction is positive
                                    }else if (tupleEnergy - lowerbound > 0) {
                                        correctionMatrix.setHigherOrder(tuple.tup, correction);
                                    } else
                                        System.err.println("Negative correction for " + tuple.tup.stringListing());
                                }

                                // Record stats
                                double delta = 2.0;
                                synchronized(sequenceBound.state){
                                    sequenceBound.state.numPartialMinimizations += result.numPartials;
                                    sequenceBound.state.totalTimePartialMin += result.timeS;
                                    sequenceBound.state.numRoundsPartialMin++;
                                    delta = sequenceBound.state.calcDelta();
                                }
                                synchronized(this){
                                    this.state.numPartialMinimizations += result.numPartials;
                                    this.state.totalTimePartialMin += result.timeS;
                                    this.state.numRoundsPartialMin++;
                                    this.progress.reportPartialMinimization(result.numPartials, delta);
                                }
                            }
                    );
                    break;
                }
                case None: {
                }

            }

            synchronized (this){
            }

        }
        loopTasks.waitForFinish();
        sequenceBound.updateBound();
        curEps = sequenceBound.state.calcDelta();
        debugPrint(String.format("Tracking Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                curEps,
                sequenceBound.state.lowerBound,
                sequenceBound.state.upperBound
        ));
        debugPrint(String.format("Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                sequenceBound.getSequenceEpsilon(),
                sequenceBound.getLowerBound(),
                sequenceBound.getUpperBound()
        ));
        debugPrint(String.format("Minimized %d nodes, %d nodes in fringe.",
                sequenceBound.finishedNodes.size(),
                sequenceBound.fringeNodes.size()));

        debugPrint(String.format("--- Minimizations --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
                this.state.numEnergiedConfs,
                this.state.totalTimeEnergy / this.state.numEnergiedConfs,
                this.state.totalTimeEnergy / this.state.numRoundsEnergy
                ));
        debugPrint(String.format("--- Expansions --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
                this.state.numExpansions,
                this.state.totalTimeExpansion / this.state.numExpansions,
                this.state.totalTimeExpansion / this.state.numRoundsExpand
        ));
        debugPrint(String.format("--- Partials --- #: %d, Avg time: %1.3e s, Avg time per round: %1.3e s",
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
        Map<BatchCorrectionMinimizer.PartialMinimizationTuple, EnergyCalculator.EnergiedParametricMolecule> confs;
        int numPartials;
        double timeS;
    }

    private void onExpansion(ExpansionResult result, SingleSequenceSHARKStarBound sequenceBound){

    }

    private void onMinimization(MinimizationResult result, SingleSequenceSHARKStarBound sequenceBound){
        // first, do some error checking
        if (result.deltaLB.compareTo(BigDecimal.ZERO) < 0)
            // print out some potentially useful information
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                result.minimizedNode.getConfSearchNode().index(context.index);
                System.err.println(String.format("Uncorrected g ub: %f, Energy %f",
                        context.partialConfUpperBoundScorer.calc(context.index, sequenceBound.seqRCs),
                        result.energy
                ));
                System.err.println(String.format("FATAL ERROR: Lower bound is decreasing with %1.9e for %s",
                        result.deltaLB,
                        Arrays.toString(result.minimizedNode.getConfSearchNode().assignments)
                        ));
            }
        if (result.deltaUB.compareTo(BigDecimal.ZERO) > 0) {
            double uncorrectedLowerBound = 0;
            double correctedLowerBound = 0;
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                result.minimizedNode.getConfSearchNode().index(context.index);
                uncorrectedLowerBound = context.partialConfLowerBoundScorer.calc(context.index, sequenceBound.seqRCs);
                correctedLowerBound = correctionMatrix.confE(result.minimizedNode.getConfSearchNode().assignments);
            }
            System.err.println(String.format("FATAL ERROR: Upper bound is increasing with %1.9e for %s",
                    result.deltaUB,
                    Arrays.toString(result.minimizedNode.getConfSearchNode().assignments)
            ));
        }

        // initialize reporting things
        double delta = 2.0;
        long fringeSize = 0;
        synchronized(sequenceBound.state) {
            // Update partition function values
            sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
            sequenceBound.state.lowerBound = sequenceBound.state.lowerBound.add(result.deltaLB, PartitionFunction.decimalPrecision);
            // Compute reporting things
            delta = sequenceBound.state.calcDelta();
            fringeSize = sequenceBound.fringeNodes.size();

            // report minimization
            sequenceBound.state.numEnergiedConfs++;
            sequenceBound.state.totalTimeEnergy+=result.timeS;
            sequenceBound.state.numRoundsEnergy++;
        }

        synchronized (this) { // don't race the main thread
            if (precomputedSequence.equals(confSpace.makeUnassignedSequence()))
                correctionMatrix.setHigherOrder(result.minimizedNode.toTuple(),
                        result.energy - minimizingEmat.confE(result.minimizedNode.getConfSearchNode().assignments));
            numConfsEnergied++;
            this.numConfsEnergiedThisLoop++;
            minList.set(result.conf.getAssignments().length - 1, minList.get(result.conf.getAssignments().length - 1) + 1);

            this.state.numEnergiedConfs++;
            this.state.totalTimeEnergy += result.timeS;
            this.state.numRoundsEnergy++;

            // report leaf
            this.progress.reportLeafNode(result.minimizedNode.getConfSearchNode().getPartialConfLowerBound(), fringeSize,delta);

            sequenceBound.addFinishedNode(result.minimizedNode);

            leafTimeAverage = result.timeS;
            debugPrint("Processed 1 leaf in " + result.timeS + " seconds.");
        }
    }

    private void onCorrection(CorrectionResult result, SingleSequenceSHARKStarBound sequenceBound){
        synchronized(sequenceBound.state){
            sequenceBound.state.upperBound = sequenceBound.state.upperBound.add(result.deltaUB, PartitionFunction.decimalPrecision);
        }
        synchronized(sequenceBound.fringeNodes){
            sequenceBound.fringeNodes.add(result.correctedNode);
        }

    }

    private void onPartialMinimization(PartialMinimizationResult result, SingleSequenceSHARKStarBound sequenceBound){

    }

    private MinimizationResult minimizeNode(MultiSequenceSHARKStarNode node, Sequence sequence, RCs seqRCs) {
        MinimizationResult result = new MinimizationResult();
        Stopwatch minTimer = new Stopwatch().start();

        // Record starting bounds
        BigDecimal startingLB = node.getLowerBound(sequence);
        BigDecimal startingUB = node.getUpperBound(sequence);


        // Get the uncorrected lower bound again for debugging purposes
        // TODO: remove the need for this extra work
        double uncorrectedLowerBound;
        try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            node.index(context.index);

            uncorrectedLowerBound = context.partialConfLowerBoundScorer.calc(context.index, seqRCs);
        }

        // Actually minimize the node
        ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.getConfSearchNode().assignments, node.getConfLowerBound(sequence));
        ConfAnalyzer.ConfAnalysis analysis = confAnalyzer.analyze(conf);

        //Submit conf for partial minimization
        energyMatrixCorrector.scheduleEnergyCorrection(analysis, conf, theBatcher);

        /* DEBUG CHECKS
        Sometimes a bug will cause the minimized energy to exceed the bounds. We must catch this to retain provable bounds.
         */
        double energy = analysis.epmol.energy;
        double newConfUpper = energy;
        double newConfLower = energy;
        // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
        double oldConfUpper = node.getConfUpperBound(sequence);
        double oldConfLower = node.getConfLowerBound(sequence);

        /* Case 1: minimized energy goes above energy upper bound
         */
        if (newConfUpper > oldConfUpper) {
            System.err.println(String.format("WARNING: Minimized energy exceeds upper bound: %s -> E: %f > UB: %f. Rejecting minimized energy.",
                    Arrays.toString(node.getConfSearchNode().assignments),
                    newConfUpper,
                    oldConfUpper
                    ));
            // Set energy to old conf upper bound
            newConfUpper = oldConfUpper;
            newConfLower = oldConfUpper;
        /* Case 2: minimized energy goes below uncorrected conf lower bound
         */
        }else if(newConfLower < uncorrectedLowerBound){
            System.err.println(String.format("WARNING: Minimized energy exceeds lower bound: %s -> E: %f < LB: %f. Rejecting minimized energy.",
                    Arrays.toString(node.getConfSearchNode().assignments),
                    newConfLower,
                    uncorrectedLowerBound
            ));
            // Set energy to old conf upper bound
            newConfUpper = uncorrectedLowerBound;
            newConfLower = uncorrectedLowerBound;
        /* Case 3: minimized energy goes below the corrected conf lower bound. This means the correction is bad

         */
        }else if(newConfLower < oldConfLower){
            System.err.println(String.format("WARNING: Bad correction: %s -> E: %f < LB: %f. Accepting minimized energy.",
                    Arrays.toString(node.getConfSearchNode().assignments),
                    newConfLower,
                    oldConfLower
            ));
            //TODO: remove correction from pool if necessary
        }

        String historyString = "";
        if (debug) {
            historyString = String.format("minimimized %s to %s from %s",
                    node.getConfSearchNode().confToString(), energy, getStackTrace());
        }
        // END DEBUG CHECKS

        node.setBoundsFromConfLowerAndUpperWithHistory(newConfLower, newConfUpper, bc.calc(newConfLower), bc.calc(newConfUpper), sequence, historyString);
        node.getConfSearchNode().setPartialConfLowerAndUpper(newConfLower, newConfUpper);
        debugPrint(String.format("Energy = %.6f, [%.6f, %.6f]",
                energy,
                oldConfLower,
                oldConfUpper
                ));
        node.markUpdated();

        // record the change in pfunc bounds
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

        return result;
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

    protected void recordCorrection(double lowerBound, double correction) {
        BigDecimal upper = bc.calc(lowerBound);
        BigDecimal corrected = bc.calc(lowerBound + correction);
        cumulativeZCorrection = cumulativeZCorrection.add(upper.subtract(corrected));
        upperReduction_PartialMin = upperReduction_PartialMin.add(upper.subtract(corrected));
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
        debugPrint("Found a leaf!");
        //rootNode.computeEpsilonErrorBounds(bound.sequence);
        nonZeroLower = true;
    }


    boolean isStable(BigDecimal stabilityThreshold, SingleSequenceSHARKStarBound seqBound) {
        return seqBound.state.numEnergiedConfs <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(seqBound.state.getUpperBound(), stabilityThreshold);
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

                child.index(context.index);
                SimpleConfSpace.Position designPos = confSpace.positions.get(nextPos);
                int nextDesignIndex = order.getNextPos(context.index, bound.seqRCs);
                SimpleConfSpace.Position nextDesignPos = null;
                if(nextDesignIndex >=0)
                    nextDesignPos = confSpace.positions.get(nextDesignIndex);
                MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild = curNode.makeOrUpdateChild(child, bound.sequence,
                        confLowerBound, confUpperBound, bc.calc(confLowerBound), bc.calc(confUpperBound), designPos, nextDesignPos, bound.seqRCs);
                MultiSequenceSHARKStarNodeChild.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, bc.calc(confLowerBound), bc.calc(confUpperBound), bound.sequence, historyString);
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
        debugPrint("Bounding "+startNode.toSeqString(bound.sequence));
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
                //System.out.println(String.format("Processed %d, %s so far. Bounds are now [%12.6e,%12.6e]",
                        //numNodes,
                        //overallLoop.getTime(2),
                        //rootNode.getLowerBound(bound.sequence),
                        //rootNode.getUpperBound(bound.sequence)));
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
                    //progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
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
                    //progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
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
                        bound.sequence, resultingLower, resultingUpper, bc.calc(resultingLower), bc.calc(resultingUpper), designPos, nextDesignPos, bound.seqRCs);
                newChild.setBoundsFromConfLowerAndUpperWithHistory(resultingLower,
                        resultingUpper, bc.calc(resultingLower), bc.calc(resultingUpper), bound.sequence, historyString);
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
