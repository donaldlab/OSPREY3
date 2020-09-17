package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.PartialConfAStarNode;
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
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.sharkstar.tools.SHARKStarEnsembleAnalyzer;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.ObjectPool;
import edu.duke.cs.osprey.tools.Stopwatch;
import javafx.util.Pair;

import java.io.BufferedWriter;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static edu.duke.cs.osprey.sharkstar.tools.MultiSequenceSHARKStarNodeStatistics.printTree;

public class MultiSequenceSHARKStarBound implements PartitionFunction {

    private final EnergyMatrixCorrector energyMatrixCorrector;
    private Sequence precomputedSequence;
    protected double targetEpsilon = 1;
    public static final boolean debug = false;
    public static final boolean suppressPrecisionWarnings = true;
    public static final boolean diveForLeaves = false;
    public static final boolean doCorrections = true;
    public static final boolean runUntilNonZero = false;
    public static final boolean skipAddingToFringe = true;
    public static double skipCutoff = 0;
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
    private final ConfIndex<PartialConfAStarNode> confIndex;
    public StaticBiggestLowerboundDifferenceOrder order;
    public final AStarPruner pruner;
    // TODO: Implement new AStarPruner for MARK*?
    protected RCs fullRCs;
    protected Parallelism parallelism;
    public ObjectPool<ScoreContext> contexts;
    private final ScorerFactory<PartialConfAStarNode> gscorerFactory;
    private final ScorerFactory<PartialConfAStarNode> rigidgscorerFactory;
    private final ScorerFactory<PartialConfAStarNode> hscorerFactory;
    private final ScorerFactory<PartialConfAStarNode> nhscorerFactory;

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

    private final Object lock = new Object();

    public double precomputedFlexComputeTime = 0; // how long did it take to compute the precomputed flex?
    public double pfuncCreationTime = 0; // how long did it take to make the singleSequence bounds?


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
        gscorerFactory = PairwiseGScorer::new;
        rigidgscorerFactory = PairwiseRigidGScorer::new;


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

        confIndex = new ConfIndex<>(rcs.getNumPos());
        rootNode = new MultiSequenceSHARKStarNode(confSpace.positions.size());
        rootNode.index(confIndex);
        double partialConfLowerbound = gscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        double partialConfUpperBound = rigidgscorerFactory.make(rigidEmat).calc(confIndex, rcs);

        rootNode.setPartialConfLowerAndUpper(partialConfLowerbound, partialConfUpperBound);

        // Initialize residue ordering
        this.order = new StaticBiggestLowerboundDifferenceOrder();
        order.setScorers(gscorerFactory.make(minimizingEmat), hscorerFactory.make(minimizingEmat));
        /* force init order */
        order.getNextPos(confIndex,fullRCs);

        // No precomputed sequence means the "precomputed" sequence is empty
        this.precomputedSequence = confSpace.makeUnassignedSequence();

        //double confLowerBound = partialConfLowerbound + hscorerFactory.make(minimizingEmat).calc(confIndex, rcs);
        //double confUpperBound = partialConfUpperBound + nhscorerFactory.make(rigidEmat).calc(confIndex, rcs);
        MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(
                partialConfLowerbound + hscorerFactory.make(minimizingEmat).calc(confIndex, rcs),
                partialConfUpperBound + nhscorerFactory.make(rigidEmat).calc(confIndex, rcs)
                );

        //rootNode.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, bc.calc(confLowerBound), bc.calc(confUpperBound), precomputedSequence, "(root initialization)");
        rootNode.setConfBounds(confBounds, this.precomputedSequence, "root initialization");
        //rootNode.setPfuncBounds(confBounds.boltzmannWeight(this.bc), this.precomputedSequence, "root initialization");
        rootNode.setErrorBound(this.bc.freeEnergy(confBounds.boltzmannWeight(this.bc).size(this.bc.mathContext)), this.precomputedSequence);

        this.contexts = new ObjectPool<>((lingored) -> {
            ScoreContext context = new ScoreContext();
            context.index = new ConfIndex<>(rcs.getNumPos());
            context.partialConfLowerBoundScorer = gscorerFactory.make(minimizingEmat);
            context.lowerBoundScorer = hscorerFactory.make(minimizingEmat);
            context.partialConfUpperBoundScorer = rigidgscorerFactory.make(rigidEmat);
            context.upperBoundScorer = nhscorerFactory.make(rigidEmat); //this is used for upper bounds, so we want it rigid
            context.ecalc = minimizingConfEcalc;

            return context;
        });
        this.theBatcher = new BatchCorrectionMinimizer(minimizingConfEcalc, minimizingEmat);

        progress = new MARKStarProgress(fullRCs.getNumPos());
        //confAnalyzer = new ConfAnalyzer(minimizingConfEcalc, minimizingEmat);
        confAnalyzer = new ConfAnalyzer(minimizingConfEcalc);
        ensembleAnalyzer = new SHARKStarEnsembleAnalyzer(minimizingEcalc, minimizingEmat);
        energyMatrixCorrector = new EnergyMatrixCorrector(this);

        setParallelism(parallelism);

        // Recording pfunc starting bounds
        //this.startLowerBound = rootNode.getLowerBound(precomputedSequence);
        //this.startUpperBound = rootNode.getUpperBound(precomputedSequence);
        this.startLowerBound = this.bc.calc(rootNode.getConfUpperBound(precomputedSequence));
        this.startUpperBound = this.bc.calc(rootNode.getConfLowerBound(precomputedSequence));
        this.minList = new ArrayList<Integer>(Collections.nCopies(rcs.getNumPos(), 0));

        // add the rootnode to the precomputed fringe
        this.precomputedFringe.add(rootNode);
    }

    private void processPrecomputedFlex(MultiSequenceSHARKStarBound precomputedFlex, SingleSequenceSHARKStarBound precomputedFlexSSBound) {
        // generate permutation matrix that will map flexible conf nodes to full conf nodes
        int[] permutationArray = genConfSpaceMapping(precomputedFlex);
        // First, change the order
        ConfIndex<PartialConfAStarNode> rootIndex = new ConfIndex<>(fullRCs.getNumPos());
        assert(this.precomputedFringe.size() == 1);
        //MultiSequenceSHARKStarNode rootNode = this.precomputedFringe.get(0);
        rootNode.index(rootIndex);
        this.order.updateForPrecomputedOrder(precomputedFlex.order, rootIndex, this.fullRCs, permutationArray);

        // Next, record information from the precomputed pfunc
        //precomputedPfunc = precomputedFlex;
        //precomputedRootNode = precomputedFlex.rootNode;
        this.precomputedSequence = precomputedFlex.confSpace.makeWildTypeSequence();
        //precomputedUpperBound = this.bc.calc(precomputedRootNode.getConfLowerBound(precomputedSequence));
        //precomputedLowerBound= this.bc.calc(precomputedRootNode.getConfUpperBound(precomputedSequence));
        precomputedUpperBound = precomputedFlexSSBound.state.getUpperBound();
        precomputedLowerBound = precomputedFlexSSBound.state.getUpperBound();

        // finally,
        //updatePrecomputedConfTree();
        updatePrecomputedFringe(precomputedFlexSSBound, permutationArray);
        mergeCorrections(precomputedFlex.correctionMatrix, genConfSpaceMapping(precomputedFlex));

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
        Stopwatch makePfuncWatch = new Stopwatch().start();
        SingleSequenceSHARKStarBound newBound = new SingleSequenceSHARKStarBound(this, seq, this);
        newBound.init(null, null, targetEpsilon);
        System.out.println("Creating new pfunc for sequence "+seq);
        System.out.println("Full RCs: "+fullRCs);
        System.out.println("Sequence RCs: "+newBound.seqRCs);
        //computeFringeForSequence(newBound, this.rootNode);

        /*
        // Compute the fringe in parallel
        Object fringeLock = new Object();
        MathTools.BigDecimalBounds fringeBounds = new MathTools.BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
        List<MultiSequenceSHARKStarNode> scoredFringe = computeFringeForSequenceParallel(newBound, this, fringeBounds, fringeLock);
        if(scoredFringe.size()==0)
            scoredFringe.add(this.rootNode);
        debugPrint(String.format("[Normal fringe # nodes, Parallel fringe # nodes] = [%d, %d]",newBound.fringeNodes.size(), scoredFringe.size()));

        //newBound.fringeNodes.addAll(scoredFringe);
        for (int i = 0; i < scoredFringe.size(); i++){
            MultiSequenceSHARKStarNode node = scoredFringe.get(i);
            if(node.getLevel() >= newBound.seqRCs.getNumPos()){
                newBound.leafQueue.add(node);
            }else{
                newBound.internalQueue.add(node);
            }
        }

        newBound.state.setBounds(fringeBounds.lower, fringeBounds.upper);
         */
        //computeFringeForSequenceParallelV2(newBound, this);
        computeRootFringe(newBound, this);
        newBound.state.updateBounds(BigDecimal.ZERO, BigDecimal.ZERO);//hack to update the state

        System.out.println(String.format("Created pfunc for %s with eps: %.6f, [%1.3e, %1.3e], %d nodes in the leaf queue, %d nodes in the internal queue",
                newBound.sequence,
                newBound.state.getDelta(),
                newBound.state.getLowerBound(),
                newBound.state.getUpperBound(),
                newBound.leafQueue.size(),
                newBound.internalQueue.size()
                ));
        if(newBound.getSequenceEpsilon() == 0)
            System.err.println(String.format("Perfectly bounded sequence for %s %s, Eps: %.9f [%1.3e, %1.3e]",
                    newBound.sequence,
                    newBound.seqRCs,
                    newBound.state.getDelta(),
                    newBound.state.getLowerBound(),
                    newBound.state.getUpperBound()
                    ));
        this.pfuncCreationTime += makePfuncWatch.stop().getTimeS();
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

    public void updatePrecomputedFringe(SingleSequenceSHARKStarBound precomputedSSBound, int[] permArray){
        final int size = this.confSpace.getNumPos();
        this.precomputedFringe.clear();
        Stream.of(precomputedSSBound.internalQueue, precomputedSSBound.leafQueue, precomputedSSBound.finishedNodes)
                .flatMap(Collection::parallelStream)
                .forEach((n) ->{
                   n.makeNodeCompatibleWithConfSpace(permArray, size);
                   this.precomputedFringe.add(n);
                });
    }

    private void updatePrecomputedNode(MultiSequenceSHARKStarNode node, int[] permutation, int size) {
        node.makeNodeCompatibleWithConfSpace(permutation, size);
        List<MultiSequenceSHARKStarNode> children = node.getAllChildren();
        if (!children.isEmpty()) {
            for (MultiSequenceSHARKStarNode child : children) {
                updatePrecomputedNode(child, permutation, size);
            }
        }else {
            precomputedFringe.add(node);
        }
    }

    private static class FringeResult{
        boolean isFringe;
        MultiSequenceSHARKStarBound msBound;
        SingleSequenceSHARKStarBound seqBound;
        List<MultiSequenceSHARKStarNode> nodes = new ArrayList<>();
        MathTools.BigDecimalBounds pfuncBounds = new MathTools.BigDecimalBounds(BigDecimal.ZERO, BigDecimal.ZERO);
    }

    private List<MultiSequenceSHARKStarNode> computeFringeForSequenceParallel(SingleSequenceSHARKStarBound singleSequencePfunc, MultiSequenceSHARKStarBound multiSequencePfunc, MathTools.BigDecimalBounds bounds, Object fringeLock) {
        final SingleSequenceSHARKStarBound seqBound = singleSequencePfunc;
        final MultiSequenceSHARKStarBound msBound = multiSequencePfunc;
        final List<MultiSequenceSHARKStarNode> scoredSeqFringe = Collections.synchronizedList(new ArrayList<>());
        final ConcurrentLinkedQueue<MultiSequenceSHARKStarNode> unscoredSeqFringe = new ConcurrentLinkedQueue<>(msBound.precomputedFringe);

        // Sometimes the precomputedFringe will be empty. It really shouldn't be, but for now just add the root Node if it is
        /*
        if(unscoredSeqFringe.isEmpty())
            unscoredSeqFringe.add(msBound.rootNode);

         */

        while(!unscoredSeqFringe.isEmpty() || loopTasks.isExpecting(1)){ // here 1 refers to the fringe computation
            final MultiSequenceSHARKStarNode node = unscoredSeqFringe.poll();
            if(node == null) // this is probably bad practice, but w/e TODO: fix this loop and make it less terrible
                continue;

            loopTasks.submitExpecting(
                ()->{
                    FringeResult result = new FringeResult();
                    try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
                        ScoreContext context = checkout.get();

                        node.index(context.index);

                        // alternative method for getting children
                        List<MultiSequenceSHARKStarNode> children = null;
                        if(node.getLevel() < seqBound.seqRCs.getNumPos()) {
                            int nextPos = msBound.order.getNextPos(context.index, seqBound.seqRCs);
                            int[] seqRcsAllowedAtPos = seqBound.seqRCs.get(nextPos);
                            children = node.getChildren(seqRcsAllowedAtPos);
                        }

                        if(children != null && !children.isEmpty()) {
                            // if there are children, just add them to queue, since we only want the fringe
                            result.isFringe = false;
                            result.nodes.addAll(children);
                        }else{
                            // If we are at a leaf, score the node
                            /* As a note, the gscores WILL NOT always be correct, because nodes will have been minimized.
                            Technically, the g LBs will be correct, but the g UBs will not. The g LBs will be dealt with
                            through the mechanism of corrections. So, just reset them
                             */
                            double gscoreLB = context.partialConfLowerBoundScorer.calc(context.index, seqBound.seqRCs);
                            double gscoreUB = context.partialConfUpperBoundScorer.calc(context.index, seqBound.seqRCs);

                            double hscoreLB = context.lowerBoundScorer.calc(context.index, seqBound.seqRCs);
                            double hscoreUB = context.upperBoundScorer.calc(context.index, seqBound.seqRCs);

                            double confLowerBound = gscoreLB + hscoreLB;
                            double confUpperBound = gscoreUB + hscoreUB;
                            double partialConfLowerBound = gscoreLB;
                            double partialConfUpperBound = gscoreUB;
                            // check if we should correct the node
                            double confCorrection = 0;
                            if(doCorrections){
                                confCorrection = msBound.correctionMatrix.confE(node.assignments);
                                if(confCorrection > gscoreLB && confCorrection < gscoreUB) {
                                    confLowerBound = confCorrection + hscoreLB;
                                    partialConfLowerBound = confCorrection;
                                }
                            }

                            MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(confLowerBound, confUpperBound);
                            /*
                            node.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound,
                                    confUpperBound,
                                    msBound.bc.calc(confLowerBound),
                                    msBound.bc.calc(confUpperBound),
                                    seqBound.sequence,
                                    historyString);
                             */
                            node.setPartialConfLowerAndUpper(partialConfLowerBound, partialConfUpperBound);
                            node.setConfBounds(confBounds,
                                    seqBound.sequence,
                                    "fringe");
                            /*
                            node.setPfuncBounds(confBounds.boltzmannWeight(msBound.bc),
                                    seqBound.sequence,
                                    "fringe");
                             */
                            MathTools.BigDecimalBounds zspaceBounds = confBounds.boltzmannWeight(msBound.bc);
                            node.setErrorBound(
                                    msBound.bc.freeEnergy(zspaceBounds.size(msBound.bc.mathContext)),
                                    seqBound.sequence);

                            String historyString = "";
                            if(debug){
                                historyString = String.format("%s: previous lower bound %f, g score %f, hscore %f, f score %f corrected score %f, from %s",
                                        node.confToString(), node.getConfLowerBound(seqBound.sequence), gscoreLB, hscoreLB, gscoreLB + hscoreLB, confCorrection, getStackTrace());
                            }

                            result.isFringe = true;
                            result.nodes.add(node);
                            result.pfuncBounds = zspaceBounds;
                            //bound.fringeNodes.add(node);
                            //scoredSeqFringe.add(node);
                        }
                    }
                    return result;
                },
                (FringeResult result)->{
                    if (result.isFringe) {
                        scoredSeqFringe.addAll(result.nodes);
                        synchronized(fringeLock){
                            bounds.upper = bounds.upper.add(result.pfuncBounds.upper, PartitionFunction.decimalPrecision);
                            bounds.lower = bounds.lower.add(result.pfuncBounds.lower, PartitionFunction.decimalPrecision);
                        }
                    }
                    else
                        unscoredSeqFringe.addAll(result.nodes);
                },
                    1 // here 1 refers to the fringe computation
            );
        }
        loopTasks.waitForFinishExpecting(1);
        return scoredSeqFringe;

    }

    private void computeRootFringe(SingleSequenceSHARKStarBound singleSequencePfunc, MultiSequenceSHARKStarBound multiSequencePfunc){
        final SingleSequenceSHARKStarBound seqBound = singleSequencePfunc;
        final MultiSequenceSHARKStarBound msBound = multiSequencePfunc;

        // set the pfunc's bounds to zero
        seqBound.state.setBoundsWithoutSideEffects(BigDecimal.ZERO, BigDecimal.ZERO);
        final MultiSequenceSHARKStarNode root = msBound.rootNode;

        try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
            ScoreContext context = checkout.get();

            root.index(context.index);

            // If we are at a leaf, score the node
                    /* As a note, the gscores WILL NOT always be correct, because nodes will have been minimized.
                    Technically, the g LBs will be correct, but the g UBs will not. The g LBs will be dealt with
                    through the mechanism of corrections. So, just reset them
                     */
            double gscoreLB = context.partialConfLowerBoundScorer.calc(context.index, seqBound.seqRCs);
            double gscoreUB = context.partialConfUpperBoundScorer.calc(context.index, seqBound.seqRCs);

            double hscoreLB = context.lowerBoundScorer.calc(context.index, seqBound.seqRCs);
            double hscoreUB = context.upperBoundScorer.calc(context.index, seqBound.seqRCs);

            double confLowerBound = gscoreLB + hscoreLB;
            double confUpperBound = gscoreUB + hscoreUB;
            double partialConfLowerBound = gscoreLB;
            double partialConfUpperBound = gscoreUB;
            // check if we should correct the node
            double confCorrection = 0;
            if(doCorrections){
                confCorrection = msBound.correctionMatrix.confE(root.assignments);
                if(confCorrection > gscoreLB && confCorrection < gscoreUB) {
                    confLowerBound = confCorrection + hscoreLB;
                    partialConfLowerBound = confCorrection;
                }
            }

            MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(confLowerBound, confUpperBound);

            root.setPartialConfLowerAndUpper(partialConfLowerBound, partialConfUpperBound);
            root.setConfBounds(confBounds,
                    seqBound.sequence,
                    "fringe");
            MathTools.BigDecimalBounds zspaceBounds = confBounds.boltzmannWeight(msBound.bc);
            root.setErrorBound(
                    msBound.bc.freeEnergy(zspaceBounds.size(msBound.bc.mathContext)),
                    seqBound.sequence);

            synchronized(lock) {
                seqBound.state.updateBoundsWithoutSideEffects(zspaceBounds.lower, zspaceBounds.upper);
                seqBound.internalQueue.add(root);
            }
        }
    }

    private void computeFringeForSequenceParallelV2(SingleSequenceSHARKStarBound singleSequencePfunc, MultiSequenceSHARKStarBound multiSequencePfunc) {
        final SingleSequenceSHARKStarBound seqBound = singleSequencePfunc;
        final MultiSequenceSHARKStarBound msBound = multiSequencePfunc;
        final ConcurrentLinkedQueue<MultiSequenceSHARKStarNode> unscoredSeqFringe = new ConcurrentLinkedQueue<>(msBound.precomputedFringe);

        // set the pfunc's bounds to zero
        seqBound.state.setBoundsWithoutSideEffects(BigDecimal.ZERO, BigDecimal.ZERO);
        // Sometimes the precomputedFringe will be empty. It really shouldn't be, but for now just add the root Node if it is
        /*
        if(unscoredSeqFringe.isEmpty())
            unscoredSeqFringe.add(msBound.rootNode);
         */

        while(!unscoredSeqFringe.isEmpty()){
            final MultiSequenceSHARKStarNode node = unscoredSeqFringe.poll();

            loopTasks.submitExpecting(
                    ()->{
                        // Compute bounds for all fringe nodes that are children of this node
                        FringeResult result = new FringeResult();
                        result.seqBound = seqBound;
                        result.msBound = msBound;
                        Queue<MultiSequenceSHARKStarNode> subtreeQueue = new LinkedList<>();
                        subtreeQueue.add(node);

                        try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
                            ScoreContext context = checkout.get();

                            while(!subtreeQueue.isEmpty()){
                                MultiSequenceSHARKStarNode curNode = subtreeQueue.poll();

                                curNode.index(context.index);

                                // Get all the existing children of this node that are compatible with this sequence
                                List<MultiSequenceSHARKStarNode> children = null;
                                if(curNode.getLevel() < seqBound.seqRCs.getNumPos()) {
                                    int nextPos = msBound.order.getNextPos(context.index, seqBound.seqRCs);
                                    int[] seqRcsAllowedAtPos = seqBound.seqRCs.get(nextPos);
                                    children = curNode.getChildren(seqRcsAllowedAtPos);
                                }

                                if(children != null && !children.isEmpty()) {
                                    // if there are children, just add them back to the queue, since we only want the fringe
                                    subtreeQueue.addAll(children);
                                }else{
                                    // If we are at a leaf, score the node
                            /* As a note, the gscores WILL NOT always be correct, because nodes will have been minimized.
                            Technically, the g LBs will be correct, but the g UBs will not. The g LBs will be dealt with
                            through the mechanism of corrections. So, just reset them
                             */
                                    double gscoreLB = context.partialConfLowerBoundScorer.calc(context.index, seqBound.seqRCs);
                                    double gscoreUB = context.partialConfUpperBoundScorer.calc(context.index, seqBound.seqRCs);

                                    double hscoreLB = context.lowerBoundScorer.calc(context.index, seqBound.seqRCs);
                                    double hscoreUB = context.upperBoundScorer.calc(context.index, seqBound.seqRCs);

                                    double confLowerBound = gscoreLB + hscoreLB;
                                    double confUpperBound = gscoreUB + hscoreUB;
                                    double partialConfLowerBound = gscoreLB;
                                    double partialConfUpperBound = gscoreUB;
                                    // check if we should correct the node
                                    double confCorrection = 0;
                                    if(doCorrections){
                                        confCorrection = msBound.correctionMatrix.confE(curNode.assignments);
                                        if(confCorrection > gscoreLB && confCorrection < gscoreUB) {
                                            confLowerBound = confCorrection + hscoreLB;
                                            partialConfLowerBound = confCorrection;
                                        }
                                    }

                                    MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(confLowerBound, confUpperBound);

                                    curNode.setPartialConfLowerAndUpper(partialConfLowerBound, partialConfUpperBound);
                                    curNode.setConfBounds(confBounds,
                                            seqBound.sequence,
                                            "fringe");
                                    MathTools.BigDecimalBounds zspaceBounds = confBounds.boltzmannWeight(msBound.bc);
                                    curNode.setErrorBound(
                                            msBound.bc.freeEnergy(zspaceBounds.size(msBound.bc.mathContext)),
                                            seqBound.sequence);

                                    result.nodes.add(curNode);
                                    result.pfuncBounds.upper = result.pfuncBounds.upper.add(zspaceBounds.upper, PartitionFunction.decimalPrecision);
                                    result.pfuncBounds.lower = result.pfuncBounds.lower.add(zspaceBounds.lower, PartitionFunction.decimalPrecision);
                                }
                            }

                        }
                        return result;
                    },
                    (FringeResult result)->{
                        // update the single sequence bound, making sure not to calc epsilon yet
                        synchronized(lock) {
                            result.seqBound.state.updateBoundsWithoutSideEffects(result.pfuncBounds.lower, result.pfuncBounds.upper);
                            // add the nodes to the appropriate single-sequence queue
                            for (int i = 0; i < result.nodes.size(); i++) {
                                MultiSequenceSHARKStarNode newNode = result.nodes.get(i);
                                if (newNode.getLevel() >= result.seqBound.seqRCs.getNumPos()) {
                                    result.seqBound.leafQueue.add(newNode);
                                } else {
                                    result.seqBound.internalQueue.add(newNode);
                                }
                            }
                        }
                    },
                    1 // here 1 refers to the fringe computation
            );
        }
        loopTasks.waitForFinishExpecting(1); //TODO: is it bad to do this??
    }

    /**
     * Generate a permutation matrix that lets us map positions from the precomputed confspace to the new confspace
     */
    public int[] genConfSpaceMapping(MultiSequenceSHARKStarBound precomputedMSBound) {
        // the permutation matrix maps confs in the precomputed flexible to the full confspace
        // Note that I think this works because Positions have equals() check the residue number
        return precomputedMSBound.confSpace.positions.stream()
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
        Stopwatch precomputeWatch = new Stopwatch().start();
        while(bound.getStatus() != Status.Estimated){
            bound.compute();
        }
        loopTasks.waitForFinish(); // we really do need this to finish before we can start on the other sequences
        this.precomputedFlexComputeTime = precomputeWatch.stop().getTimeS();
        if(doCorrections)
            addFullMinimizationsToCorrectionMatrix(precompFlex, bound);
        precompFlex.printEnsembleAnalysis();
        processPrecomputedFlex(precompFlex, bound);
        // Check to make sure bounds are the same as the queue bounds
        debugPrint(String.format("State eps: %.9f, [%1.9e, %1.9e]",
                bound.state.getDelta(),
                bound.state.getLowerBound(),
                bound.state.getUpperBound()
                ));
        debugPrint(String.format("Queue eps: %.9f, [%1.9e, %1.9e]",
                bound.getEpsFromQueues(),
                bound.getLowerFromQueues(),
                bound.getUpperFromQueues()
        ));
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
       SingleSequenceSHARKStarBound sequenceBound;
    }

    /**
     * Try to apply corrections to node.
     */
    private CorrectionResult correctNodeOrFalse(MultiSequenceSHARKStarNode node, SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarBound msBound) {
        CorrectionResult result = new CorrectionResult();
        result.sequenceBound = bound;
        double confCorrection = msBound.correctionMatrix.confE(node.assignments);
        double oldg = node.getPartialConfLowerBound();
        double correctionDiff = confCorrection - oldg;

        // Make sure to check that we do not overcorrect
        boolean overcorrecting = false;
        double oldGUpper = node.getPartialConfUpperBound();
        if(confCorrection > oldGUpper) {
            System.err.printf("Attempted overcorrection of %s for seq %s: [%.9f, %.9f] -> + %.9f -> [%.9f, %.9f]",
                    Arrays.toString(node.assignments),
                    bound.sequence,
                    oldg,
                    oldGUpper,
                    correctionDiff,
                    confCorrection,
                    oldGUpper
                    );
            overcorrecting = true;
        }

        if ( correctionDiff > 1e-5 && !overcorrecting && doCorrections) {
            result.didCorrect = true;
            result.correctionSize = correctionDiff;

            double oldConfLowerBound = node.getConfLowerBound(bound.sequence);
            //BigDecimal oldZUpperBound = node.getUpperBound(bound.sequence);
            BigDecimal oldZUpperBound = msBound.bc.calc(oldConfLowerBound);

            // update the node gscore
            node.setPartialConfLowerAndUpper(confCorrection, node.getPartialConfUpperBound());
            recordCorrection(oldg, correctionDiff);
            String historyString = String.format("%s: correction from %f to %f, from ",
                    node.confToString(), oldg, confCorrection);

            // update the node total scores
            MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(
                    oldConfLowerBound + correctionDiff,
                    node.getConfUpperBound(bound.sequence));
            node.setConfBounds(confBounds,
                    bound.sequence,
                    "correction");
            /*
            node.setPfuncBounds(confBounds.boltzmannWeight(msBound.bc),
                    bound.sequence,
                    "correction");
             */
            node.setErrorBound(
                    msBound.bc.freeEnergy(confBounds.boltzmannWeight(msBound.bc).size(msBound.bc.mathContext)),
                    bound.sequence);
            //node.setBoundsFromConfLowerAndUpperWithHistory(oldConfLowerBound + correctionDiff, node.getConfUpperBound(bound.sequence), bc.calc(oldConfLowerBound+correctionDiff), msBound.bc.calc(node.getConfUpperBound(bound.sequence)), bound.sequence, historyString);
            debugPrint("Correcting " + node.toSeqString(bound.sequence) +" correction ="+(correctionDiff) );

            //BigDecimal newUpperBound = node.getUpperBound(bound.sequence);
            BigDecimal newUpperBound = msBound.bc.calc(node.getConfLowerBound(bound.sequence));
            BigDecimal diffUB = newUpperBound.subtract(oldZUpperBound, PartitionFunction.decimalPrecision);
            if(diffUB.compareTo(BigDecimal.ZERO) > 0){
                throw new RuntimeException();
            }
            //if(node.getUpperBound(bound.sequence).compareTo(node.getLowerBound(bound.sequence)) < 0){
            if(node.getConfLowerBound(bound.sequence) > node.getConfUpperBound(bound.sequence)){
                System.err.println(String.format("Insane bounds: LB (%.3f) > LB (%.3f)",
                        node.getConfLowerBound(bound.sequence),
                        node.getConfUpperBound(bound.sequence)
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
        SingleSequenceSHARKStarBound sequenceBound;
        MultiSequenceSHARKStarBound msBound;
    }

    public void computeForSequenceParallel(int maxNumConfs, SingleSequenceSHARKStarBound sequenceBound){
        debugPrint("Tightening bound for "+sequenceBound.sequence + " " + this.cachePattern);
        Stopwatch computeWatch = new Stopwatch().start();

        double lastEps = sequenceBound.state.getDelta();

        /*
        if(lastEps == 0)
            System.err.println("Epsilon is ZERO??! And we are still tightening the bound!?");
         */

        long previousConfCount = sequenceBound.state.workDone();

        // Run to populate the queues a little bit -- not clear that this is a good thing, since it's not parallelized, but hey
        if (!sequenceBound.nonZeroLower() && runUntilNonZero){
            runUntilNonZero(sequenceBound, this);
            sequenceBound.updateStateFromQueues();
        }
        Step step = Step.None;

        /*
        Continue looping if 1) there are nodes in the fringe queue or
        2) we are running tasks that will result in nodes being added to the queue
         */
        this.numConfsEnergiedThisLoop = 0;
        int numRoundsWorking = 0;

        computation: while( (!sequenceBound.internalQueue.isEmpty() || !sequenceBound.leafQueue.isEmpty() || loopTasks.isExpecting(0))){ // here, 0 refers to the pfunc computation
            synchronized(lock) {
                step = Step.None;
                // Record starting delta
                double newEps = sequenceBound.state.getDelta();
                if(lastEps - newEps < -1e13)
                    System.err.println(String.format("ERROR: Epsilon increased from %f to %f, difference of %1.3e",
                            lastEps, newEps, lastEps-newEps));
                lastEps = newEps;

                debugPrint(String.format("Epsilon: %.9f, Bounds:[%1.3e, %1.3e]",
                        lastEps,
                        sequenceBound.state.getLowerBound(),
                        sequenceBound.state.getUpperBound()
                ));

                // Early termination
                if (/*lastEps < targetEpsilon ||*/ sequenceBound.getStatus() == Status.Estimated ||
                        //sequenceBound.state.workDone() - previousConfCount >= maxNumConfs ||
                        numRoundsWorking >= maxNumConfs ||
                        !isStable(stabilityThreshold, sequenceBound)
                ){
                    if(sequenceBound.state.workDone() - previousConfCount >= maxNumConfs)
                        debugPrint("Exiting loop because of work");

                    if(!isStable(stabilityThreshold, sequenceBound)) {
                        debugPrint("Exiting loop due to stablity, thresh: " + stabilityThreshold);
                        sequenceBound.setStatus(Status.Unstable);
                    }
                    break;
                }
            }


            // Make variables to pass to the child threads
            List<MultiSequenceSHARKStarNode> internalNodes = new ArrayList<>();
            MultiSequenceSHARKStarNode loosestLeaf = null;
            BatchCorrectionMinimizer.Batch batch = null;

            // Make stopwatches
            Stopwatch loopWatch = new Stopwatch();
            loopWatch.start();
            Stopwatch internalTime = new Stopwatch();

            double expansionsPerEnergy = (sequenceBound.state.totalTimeEnergy / sequenceBound.state.numEnergiedConfs) /
                    (sequenceBound.state.totalTimeExpansion / sequenceBound.state.numExpansions);
            int maxNumNodes = Math.min((int) (expansionsPerEnergy * 1e-1), 1000);
            if(maxNumNodes < 1)
                maxNumNodes=1;
            int numNodes = 0;

            synchronized(lock) {
                //Figure out what step to make and get nodes
                boolean computePartials = theBatcher.canProcess();
                // If we have partial minimizations to do, do them with highest priority
                if (computePartials && doCorrections) {
                    batch = theBatcher.getBatch();
                    step = Step.Partial;
                // If we have nodes in either the leaf or internal queues, process them
                } else if (sequenceBound.leafQueue.size() > 0 || sequenceBound.internalQueue.size() > 0) {
                    //Record the largest error in the leaf queue and internal queues
                    double internalErrorWeightFactor = 1e3;
                    double freeEnergyLeafError = Double.POSITIVE_INFINITY;
                    double freeEnergyInternalError = Double.POSITIVE_INFINITY;
                    if (!sequenceBound.leafQueue.isEmpty())
                        freeEnergyLeafError = sequenceBound.leafQueue.peek().getErrorBound(sequenceBound.sequence);
                    if (!sequenceBound.internalQueue.isEmpty())
                        freeEnergyInternalError = (sequenceBound.internalQueue.peek().getErrorBound(sequenceBound.sequence)
                        - BoltzmannCalculator.constRT*Math.log(internalErrorWeightFactor));

                    // If the leaf error is larger then (less than since we are in free energy space) the internal error,
                    // then grab the leaf and do an energy step
                    if (freeEnergyLeafError < Double.POSITIVE_INFINITY && freeEnergyLeafError < freeEnergyInternalError) {
                            loosestLeaf = sequenceBound.leafQueue.poll();
                            step = Step.Energy;
                    // If the internal error is larger then (less than since we are in free energy space) the leaf error,
                    //then grab the internal node and do an expansion step
                    } else if (freeEnergyInternalError < Double.POSITIVE_INFINITY && freeEnergyInternalError <= freeEnergyLeafError) {
                        double nodeError = freeEnergyInternalError;
                        while(nodeError < Double.POSITIVE_INFINITY &&
                                nodeError <= freeEnergyLeafError &&
                                internalNodes.size() < maxNumNodes
                        ){
                            internalNodes.add(sequenceBound.internalQueue.poll());
                            if(!sequenceBound.internalQueue.isEmpty())
                                nodeError = (sequenceBound.internalQueue.peek().getErrorBound(sequenceBound.sequence)
                                        - BoltzmannCalculator.constRT*Math.log(internalErrorWeightFactor));
                            else
                                nodeError = Double.POSITIVE_INFINITY;
                        }

                        step = Step.ExpandInBatches;
                    // Otherwise, do nothing and print an error
                    } else {
                        step = Step.None;
                        System.err.println(String.format("No step even though we have nodes: # internals: %d, # leaves: %d, internal error: %.3f, leaf error %.3f",
                                sequenceBound.internalQueue.size(),
                                sequenceBound.leafQueue.size(),
                                freeEnergyInternalError,
                                freeEnergyLeafError
                                ));
                    }
                }else{
                    step = Step.None;
                }
            }


            // Take steps
            switch(step) {
                case Energy:{
                    // Make the variables final so that they can be somewhat functional
                    MultiSequenceSHARKStarNode toMinimize = loosestLeaf;
                    SingleSequenceSHARKStarBound singleSequencePfunc = sequenceBound;
                    MultiSequenceSHARKStarBound multiSequencePfunc = this;

                    debugPrint(String.format("Minimizing node with lower bound %f and pfunc error %1.3e",
                            loosestLeaf.getConfLowerBound(singleSequencePfunc.sequence), //TODO: fix the null-pointer exception thrown here
                            loosestLeaf.getErrorBound(singleSequencePfunc.sequence)
                            ));

                    loopTasks.submit( () -> {
                        CorrectionResult correctionResult = correctNodeOrFalse(toMinimize, singleSequencePfunc, multiSequencePfunc);
                        if(correctionResult.didCorrect)
                            return correctionResult;
                        else
                            return minimizeNode(toMinimize, singleSequencePfunc, multiSequencePfunc);
                            },
                            (result) -> {
                        if (result.getClass() == CorrectionResult.class)
                            onCorrection((CorrectionResult) result);
                        else
                            onMinimization((MinimizationResult) result);
                            }
                    );
                    numRoundsWorking++;
                    break;
                }
                case ExpandInBatches: {
                    final List<MultiSequenceSHARKStarNode> toExpand = internalNodes;
                    final SingleSequenceSHARKStarBound singleSequencePfunc = sequenceBound;
                    final MultiSequenceSHARKStarBound multiSequencePfunc = this;
                    numNodes = toExpand.size();
                    debugPrint("Processing " + numNodes + " internal nodes...");

                    loopTasks.submitExpecting(
                            () -> {
                                ExpansionResult result = new ExpansionResult();
                                result.sequenceBound = singleSequencePfunc;
                                result.msBound = multiSequencePfunc;
                                internalTime.reset();
                                internalTime.start();

                                BigDecimal startLB = BigDecimal.ZERO;
                                BigDecimal startUB = BigDecimal.ZERO;

                                result.newNodes = Collections.synchronizedList(new ArrayList<>());
                                for (MultiSequenceSHARKStarNode internalNode : toExpand) {

                                    startLB = startLB.add(
                                            multiSequencePfunc.bc.calc(internalNode.getConfUpperBound(singleSequencePfunc.sequence)),
                                            PartitionFunction.decimalPrecision);
                                    startUB = startUB.add(multiSequencePfunc.bc.calc(internalNode.getConfLowerBound(singleSequencePfunc.sequence)),
                                            PartitionFunction.decimalPrecision);

                                    if (diveForLeaves && !MathTools.isGreaterThan(
                                            multiSequencePfunc.bc.calc(internalNode.getConfUpperBound(singleSequencePfunc.sequence)),
                                            BigDecimal.ONE) &&
                                            (singleSequencePfunc.state.getUpperBound().compareTo(BigDecimal.ONE) < 0 ||
                                            MathTools.isGreaterThan(
                                                    MathTools.bigDivide(
                                                            multiSequencePfunc.bc.calc(internalNode.getConfLowerBound(singleSequencePfunc.sequence)),
                                                            singleSequencePfunc.state.getUpperBound(),
                                                            PartitionFunction.decimalPrecision),
                                                    new BigDecimal(1 - multiSequencePfunc.targetEpsilon)))
                                    ) {
                                        boundLowestBoundConfUnderNode(singleSequencePfunc, multiSequencePfunc, internalNode, result.newNodes);
                                    } else {

                                        processPartialConfNode(singleSequencePfunc, multiSequencePfunc, result.newNodes, internalNode);
                                    }
                                }

                                BigDecimal endLB = result.newNodes.stream()
                                        .map( n -> multiSequencePfunc.bc.calc(n.getConfUpperBound(singleSequencePfunc.sequence)))
                                        .reduce(BigDecimal.ZERO, (a,b) -> a.add(b, PartitionFunction.decimalPrecision));
                                BigDecimal endUB = result.newNodes.stream()
                                        .map( n -> multiSequencePfunc.bc.calc(n.getConfLowerBound(singleSequencePfunc.sequence)))
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
                                            toExpand.get(0).dumpHistory(singleSequencePfunc.sequence);
                                            toExpand.get(0).getAllChildren().forEach((n)-> n.dumpHistory(singleSequencePfunc.sequence));
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
                            this::onExpansion,
                            1// here, 0 refers to the pfunc computation
                    );
                    numRoundsWorking++;
                    break;
                }
                case Partial: {
                    final SingleSequenceSHARKStarBound singleSequencePfunc = sequenceBound;
                    final MultiSequenceSHARKStarBound multiSequencePfunc = this;
                    final BatchCorrectionMinimizer.Batch theFinalBatch = batch;
                    debugPrint("Computing partial mins");
                    loopTasks.submit(
                            () -> {
                                PartialMinimizationResult result = new PartialMinimizationResult();
                                result.sequenceBound = singleSequencePfunc;
                                result.msBound = multiSequencePfunc;
                                Stopwatch partialMinTime = new Stopwatch().start();

                                // calculate all the fragment energies
                                Map<BatchCorrectionMinimizer.PartialMinimizationTuple, EnergyCalculator.EnergiedParametricMolecule> confs = new HashMap<>();
                                for (int i =0; i< theFinalBatch.fragments.size(); i++) {
                                    BatchCorrectionMinimizer.PartialMinimizationTuple frag = theFinalBatch.fragments.get(i);
                                    double energy;

                                    // are there any RCs are from two different backbone states that can't connect?
                                    if (result.msBound.theBatcher.isParametricallyIncompatible(frag.tup)) {

                                        // yup, give this frag an infinite energy so we never choose it
                                        energy = Double.POSITIVE_INFINITY;

                                    } else {

                                        // nope, calculate the usual fragment energy
                                        confs.put(frag, result.msBound.theBatcher.confEcalc.calcEnergy(frag.tup));
                                    }
                                }
                                partialMinTime.stop();

                                result.confs = confs;
                                result.numPartials = confs.size();
                                result.timeS = partialMinTime.getTimeS();
                                return result;
                            },
                            this::onPartialMinimization
                    );
                    break;
                }
                case None: {
                    // This case happens when we have no good nodes to process
                    System.err.println(String.format("None step for state: %s %s, Eps: %.6f [%1.9e, %1.9e]",
                            sequenceBound.sequence,
                            this.cachePattern,
                            sequenceBound.state.getDelta(),
                            sequenceBound.state.getLowerBound(),
                            sequenceBound.state.getUpperBound()
                            ));
                    break computation;
                }

            }

        }
        //loopTasks.waitForFinish();
        //sequenceBound.updateBound();
        double timeElapsed = computeWatch.stop().getTimeS();

        if(this.state.secondsPerSeq.containsKey(sequenceBound.sequence))
            this.state.secondsPerSeq.replace(sequenceBound.sequence, this.state.secondsPerSeq.get(sequenceBound.sequence) + timeElapsed);
        else
            this.state.secondsPerSeq.put(sequenceBound.sequence, timeElapsed);

        lastEps = sequenceBound.state.getDelta();
        debugPrint(String.format("Tracking Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                lastEps,
                sequenceBound.state.getLowerBound(),
                sequenceBound.state.getUpperBound()
        ));
        debugPrint(String.format("Queue Epsilon: %.9f, Bounds:[%1.9e, %1.9e]",
                sequenceBound.getEpsFromQueues(),
                sequenceBound.getLowerFromQueues(),
                sequenceBound.getUpperFromQueues()
        ));
        debugPrint(String.format("Minimized %d nodes, %d nodes in fringe.",
                sequenceBound.finishedNodes.size(),
                sequenceBound.internalQueue.size() + sequenceBound.leafQueue.size()));

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
        SingleSequenceSHARKStarBound sequenceBound; // Pass this through to try to remove synchronization requirement
        MultiSequenceSHARKStarBound msBound; // Pass this through to try to remove synchronization requirement
    }

    private static class PartialMinimizationResult{
        Map<BatchCorrectionMinimizer.PartialMinimizationTuple, EnergyCalculator.EnergiedParametricMolecule> confs;
        int numPartials;
        double timeS;
        SingleSequenceSHARKStarBound sequenceBound;
        MultiSequenceSHARKStarBound msBound;
    }

    private void onExpansion(ExpansionResult result){
        double delta = 2.0;
        long fringeSize = 0;
        synchronized(result.msBound.lock) {
            // Update partition function values
            result.sequenceBound.state.updateBounds(result.deltaLB, result.deltaUB);
            result.sequenceBound.state.numExpansions += result.numExpanded;
            result.sequenceBound.state.totalTimeExpansion += result.timeS;
            result.sequenceBound.state.numRoundsExpand++;
            delta = result.sequenceBound.state.getDelta();
            fringeSize = result.sequenceBound.internalQueue.size() + result.sequenceBound.leafQueue.size();

            // update the tracking variables
            result.msBound.state.numExpansions += result.numExpanded;
            result.msBound.state.totalTimeExpansion += result.timeS;
            result.msBound.state.numRoundsExpand++;
            debugPrint("Got " + result.newNodes.size() +" internal nodes.");
            result.msBound.internalTimeSum = result.timeS;
            result.msBound.internalTimeAverage = result.msBound.internalTimeSum / Math.max(1, result.numExpanded);
            debugPrint("Internal node time :" + result.msBound.internalTimeSum + ", average " + result.msBound.internalTimeAverage);
            result.msBound.numInternalNodesProcessed += result.numExpanded;
            // report the nodes
            for (int i = 0; i< result.newNodes.size(); i++){
                MultiSequenceSHARKStarNode node = result.newNodes.get(i);
                result.msBound.progress.reportInternalNode(node.getLevel(),
                        node.getGScore(),
                        node.getConfLowerBound(result.sequenceBound.sequence) - node.getGScore(),
                        fringeSize,
                        result.newNodes.size(),
                        delta);
            }
            // add the nodes back into the queue
            for (int i = 0; i < result.newNodes.size(); i++){
                MultiSequenceSHARKStarNode node = result.newNodes.get(i);
                if (
                        node.getErrorBound(result.sequenceBound.sequence) < skipCutoff ||
                        result.sequenceBound.sequence.equals(this.precomputedSequence) ||
                        !skipAddingToFringe
                        ) {
                    if (node.getLevel() >= result.sequenceBound.seqRCs.getNumPos()) {
                        result.sequenceBound.leafQueue.add(node);
                    } else {
                        result.sequenceBound.internalQueue.add(node);
                    }
                }else{
                    debugPrint(String.format("Skipping fringe addition of node with no error"));
                }
            }
        }

    }

    private void onMinimization(MinimizationResult result){
        // Try to make a batch of partial minimizations TODO: should this go here or elsewhere?
        while(result.msBound.theBatcher.canBatch())
            result.msBound.theBatcher.makeBatch();
        // do some error checking
        if (result.deltaLB.compareTo(BigDecimal.ZERO) < 0)
            // print out some potentially useful information
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                result.minimizedNode.index(context.index);
                System.err.println(String.format("Uncorrected g ub: %f, Energy %f",
                        context.partialConfUpperBoundScorer.calc(context.index, result.sequenceBound.seqRCs),
                        result.energy
                ));
                System.err.println(String.format("FATAL ERROR: Lower bound is decreasing with ΔZ %1.9e for %s",
                        result.deltaLB,
                        Arrays.toString(result.minimizedNode.assignments)
                ));
            }
        if (result.deltaUB.compareTo(BigDecimal.ZERO) > 0) {
            double uncorrectedLowerBound = 0;
            double correctedLowerBound = 0;
            try (ObjectPool.Checkout<ScoreContext> checkout = contexts.autoCheckout()) {
                ScoreContext context = checkout.get();
                result.minimizedNode.index(context.index);
                uncorrectedLowerBound = context.partialConfLowerBoundScorer.calc(context.index, result.sequenceBound.seqRCs);
                correctedLowerBound = correctionMatrix.confE(result.minimizedNode.assignments);
            }
            System.err.println(String.format("FATAL ERROR: Upper bound is increasing with ΔZ %1.9e for %s",
                    result.deltaUB,
                    Arrays.toString(result.minimizedNode.assignments)
            ));
        }

        // initialize reporting things
        double delta = 2.0;
        long fringeSize = 0;
        synchronized(result.msBound.lock) {
            // Update partition function values
            result.sequenceBound.state.updateBounds(result.deltaLB, result.deltaUB);
            // Compute reporting things
            delta = result.sequenceBound.state.getDelta();
            fringeSize = result.sequenceBound.internalQueue.size() + result.sequenceBound.leafQueue.size();

            // report minimization
            result.sequenceBound.state.numEnergiedConfs++;
            result.sequenceBound.state.totalTimeEnergy+=result.timeS;
            result.sequenceBound.state.numRoundsEnergy++;

            //update tracking variables
            result.msBound.numConfsEnergied++;
            result.msBound.numConfsEnergiedThisLoop++;
            result.msBound.minList.set(result.conf.getAssignments().length - 1, result.msBound.minList.get(result.conf.getAssignments().length - 1) + 1);

            result.msBound.state.numEnergiedConfs++;
            result.msBound.state.totalTimeEnergy += result.timeS;
            result.msBound.state.numRoundsEnergy++;

            // report leaf
            result.msBound.progress.reportLeafNode(result.minimizedNode.getPartialConfLowerBound(), fringeSize,delta);

            result.sequenceBound.addFinishedNode(result.minimizedNode);

            result.msBound.leafTimeAverage = result.timeS;
            debugPrint(String.format("Processed %s: [%.3f, ??] -> %.3f in %.3f seconds.",
                    result.minimizedNode.confToString(),
                    result.conf.getScore(),
                    result.energy,
                    result.timeS
                    ));
        }
    }

    private void onCorrection(CorrectionResult result){
        synchronized(lock){
            result.sequenceBound.state.updateBounds(BigDecimal.ZERO, result.deltaUB);

            result.sequenceBound.leafQueue.add(result.correctedNode);
        }

    }

    private void onPartialMinimization(PartialMinimizationResult result){
        // update the energy matrix
        for(BatchCorrectionMinimizer.PartialMinimizationTuple tuple : result.confs.keySet()) {
            double lowerbound = result.msBound.theBatcher.minimizingEnergyMatrix.getInternalEnergy(tuple.tup);
            double uncorrectedParentConfLB = result.msBound.theBatcher.minimizingEnergyMatrix.getInternalEnergy(tuple.parentConf);
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
                result.msBound.correctionMatrix.setHigherOrder(tuple.tup, correction);
            } else
                System.err.println("Negative correction for " + tuple.tup.stringListing());
        }

        // Record stats
        double delta = 2.0;
        synchronized(lock){
            result.sequenceBound.state.numPartialMinimizations += result.numPartials;
            result.sequenceBound.state.totalTimePartialMin += result.timeS;
            result.sequenceBound.state.numRoundsPartialMin++;
            delta = result.sequenceBound.state.getDelta();

            // update tracking vars
            result.msBound.state.numPartialMinimizations += result.numPartials;
            result.msBound.state.totalTimePartialMin += result.timeS;
            result.msBound.state.numRoundsPartialMin++;
            result.msBound.progress.reportPartialMinimization(result.numPartials, delta);
        }
    }

    private MinimizationResult minimizeNode(MultiSequenceSHARKStarNode node, SingleSequenceSHARKStarBound sequenceBound, MultiSequenceSHARKStarBound msBound) {
        MinimizationResult result = new MinimizationResult();
        result.sequenceBound = sequenceBound;
        result.msBound = msBound;
        Stopwatch minTimer = new Stopwatch().start();

        // Record starting bounds
        BigDecimal startingLB = msBound.bc.calc(node.getConfUpperBound(result.sequenceBound.sequence));
        BigDecimal startingUB = msBound.bc.calc(node.getConfLowerBound(result.sequenceBound.sequence));


        // Get the uncorrected lower bound again for debugging purposes
        // TODO: remove the need for this extra work
        double uncorrectedLowerBound;
        try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            node.index(context.index);

            uncorrectedLowerBound = context.partialConfLowerBoundScorer.calc(context.index, result.sequenceBound.seqRCs);
        }

        // Actually minimize the node
        ConfSearch.ScoredConf conf = new ConfSearch.ScoredConf(node.assignments, node.getConfLowerBound(result.sequenceBound.sequence));
        ConfAnalyzer.ConfAnalysis analysis = msBound.confAnalyzer.analyze(conf);

        //Submit conf for partial minimization
        msBound.energyMatrixCorrector.scheduleEnergyCorrection(analysis, conf, theBatcher);

        /* DEBUG CHECKS
        Sometimes a bug will cause the minimized energy to exceed the bounds. We must catch this to retain provable bounds.
         */
        double energy = analysis.epmol.energy;
        double newConfUpper = energy;
        double newConfLower = energy;
        // Record pre-minimization bounds so we can parse out how much minimization helped for upper and lower bounds
        double oldConfUpper = node.getConfUpperBound(result.sequenceBound.sequence);
        double oldConfLower = node.getConfLowerBound(result.sequenceBound.sequence);

        /* Case 1: minimized energy goes above energy upper bound
         */
        if (newConfUpper > oldConfUpper) {
            System.err.println(String.format("WARNING: Minimized energy exceeds upper bound: %s -> E: %f > UB: %f. Rejecting minimized energy.",
                    Arrays.toString(node.assignments),
                    newConfUpper,
                    oldConfUpper
                    ));
            // Set energy to old conf upper bound
            newConfUpper = oldConfUpper;
            newConfLower = oldConfUpper;
        /* Case 2: minimized energy goes below the corrected conf lower bound, but not the uncorrected. This means the correction is bad

         */
        }else if(newConfLower < oldConfLower && newConfLower > uncorrectedLowerBound){
            System.err.println(String.format("WARNING: Bad correction: %s -> E: %f < LB: %f. Accepting minimized energy.",
                    Arrays.toString(node.assignments),
                    newConfLower,
                    oldConfLower
            ));
            //TODO: remove correction from pool if necessary
        /* Case 3: minimized energy goes below uncorrected conf lower bound
         */
        }else if(newConfLower < uncorrectedLowerBound){
            System.err.println(String.format("WARNING: Minimized energy exceeds lower bound: %s -> E: %f < LB: %f. Rejecting minimized energy.",
                    Arrays.toString(node.assignments),
                    newConfLower,
                    uncorrectedLowerBound
            ));
            // Set energy to old conf lower bound
            newConfUpper = uncorrectedLowerBound;
            newConfLower = uncorrectedLowerBound;
        }

        String historyString = "";
        if (msBound.debug) {
            historyString = String.format("minimimized %s to %s from %s",
                    node.confToString(), energy, getStackTrace());
        }
        // END DEBUG CHECKS

        MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(newConfLower, newConfUpper);
        node.setConfBounds(confBounds,
                sequenceBound.sequence,
                "minimization");
        /*
        node.setPfuncBounds(confBounds.boltzmannWeight(msBound.bc),
                sequenceBound.sequence,
                "minimization");
         */
        node.setErrorBound(
                msBound.bc.freeEnergy(confBounds.boltzmannWeight(msBound.bc).size(msBound.bc.mathContext)),
                sequenceBound.sequence);
        //node.setBoundsFromConfLowerAndUpperWithHistory(newConfLower, newConfUpper, msBound.bc.calc(newConfLower), msBound.bc.calc(newConfUpper), result.sequenceBound.sequence, historyString);
        node.setPartialConfLowerAndUpper(newConfLower, newConfUpper);
        debugPrint(String.format("Energy = %.6f, [%.6f, %.6f]",
                energy,
                oldConfLower,
                oldConfUpper
                ));

        // record the change in pfunc bounds
        BigDecimal endLB = msBound.bc.calc(node.getConfUpperBound(result.sequenceBound.sequence));
        BigDecimal endUB = msBound.bc.calc(node.getConfLowerBound(result.sequenceBound.sequence));

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
    private void runUntilNonZero(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarBound msBound) {
        System.out.println("Running until leaf is found...");
        List<MultiSequenceSHARKStarNode> newNodes = new ArrayList<>();
        if(bound.internalQueue.isEmpty())
            return;
        boundLowestBoundConfUnderNode(bound, msBound, bound.internalQueue.poll(),newNodes);
        for(MultiSequenceSHARKStarNode newNode : newNodes) {
            if(!newNode.isMinimized(bound.sequence))
                if(newNode.getLevel() >= bound.seqRCs.getNumPos()){
                    bound.leafQueue.add(newNode);
                }else{
                    bound.internalQueue.add(newNode);
                }
            else {
                bound.addFinishedNode(newNode);
            }
        }
        newNodes.clear();
        debugPrint("Found a leaf!");
        msBound.nonZeroLower = true;
    }


    boolean isStable(BigDecimal stabilityThreshold, SingleSequenceSHARKStarBound seqBound) {
        return seqBound.state.numEnergiedConfs <= 0 || stabilityThreshold == null
                || MathTools.isGreaterThanOrEqual(seqBound.state.getUpperBound(), stabilityThreshold);
    }

    private MultiSequenceSHARKStarNode drillDown(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarBound msBound,
                                                 List<MultiSequenceSHARKStarNode> newNodes,
                                                 MultiSequenceSHARKStarNode curNode) {
        PriorityQueue<MultiSequenceSHARKStarNode> queue = bound.fringeNodes;
        RCs RCs = bound.seqRCs;
        try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
            ScoreContext context = checkout.get();
            curNode.index(context.index);
            // which pos to expand next?
            int nextPos = msBound.order.getNextPos(context.index, RCs);
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
                MultiSequenceSHARKStarNode child;
                if(curNode.hasExistingChild(nextRc)){
                    child = curNode.getExistingChild(nextRc);
                }else{
                    child = curNode.assign(nextPos, nextRc);
                    curNode.addChild(child, bound.sequence);
                }
                double confLowerBound = Double.POSITIVE_INFINITY;
                double confUpperBound = Double.NEGATIVE_INFINITY;
                String historyString = "Error!";
                curNode.index(context.index);

                // score the child node differentially against the parent node
                if (child.getLevel() < RCs.getNumPos()) {

                    double confCorrection = msBound.correctionMatrix.confE(child.assignments);
                    double diff = confCorrection;
                    double rigiddiff = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double hdiff = context.lowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    double maxhdiff = context.upperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    //Correct for incorrect gscore.
                    //rigiddiff = rigiddiff - node.partialConfLowerbound + node.partialConfUpperBound;
                    child.setPartialConfLowerAndUpper(diff, rigiddiff);

                    confLowerBound = child.getPartialConfLowerBound() + hdiff;
                    confUpperBound = rigiddiff + maxhdiff;
                    if (diff < confCorrection) {
                        recordCorrection(confLowerBound, confCorrection - diff);
                        confLowerBound = confCorrection + hdiff;
                    }

                    historyString = String.format("%s: previous lower bound (none), g score %f, hscore %f, f score %f corrected score %f, from %s",
                            child.confToString(), diff, hdiff, diff+hdiff, confCorrection, getStackTrace());
                    msBound.progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
                }
                if (child.getLevel() == RCs.getNumPos()) {
                    double confRigid = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                    //confRigid = confRigid - node.partialConfLowerbound + node.partialConfUpperBound;

                    double confCorrection = msBound.correctionMatrix.confE(child.assignments);
                    double lowerbound = msBound.minimizingEmat.confE(child.assignments);
                    if (lowerbound < confCorrection) {
                        recordCorrection(lowerbound, confCorrection - lowerbound);
                    }
                    checkBounds(confCorrection, confRigid);
                    historyString = String.format("%s: pairwise lower bound: %f, previous lower bound (none), correctedbound %f, from %s",
                            child.confToString(), lowerbound, confCorrection, getStackTrace());
                    confLowerBound = confCorrection;
                    child.setPartialConfLowerAndUpper(confCorrection, confRigid);
                    confUpperBound = confRigid;
                    child.setPartialConfLowerAndUpper(confLowerBound, confUpperBound);
                    msBound.numConfsScored++;
                    msBound.progress.reportLeafNode(child.getPartialConfLowerBound(), queue.size(), bound.getSequenceEpsilon());
                }
                partialTime.stop();

                child.index(context.index);
                SimpleConfSpace.Position designPos = msBound.confSpace.positions.get(nextPos);
                int nextDesignIndex = msBound.order.getNextPos(context.index, bound.seqRCs);
                SimpleConfSpace.Position nextDesignPos = null;
                if(nextDesignIndex >=0)
                    nextDesignPos = msBound.confSpace.positions.get(nextDesignIndex);
                //MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild = curNode.makeOrUpdateChild(child, bound.sequence,
                        //confLowerBound, confUpperBound, msBound.bc.calc(confLowerBound), msBound.bc.calc(confUpperBound), designPos, nextDesignPos, bound.seqRCs);
                MultiSequenceSHARKStarNode MultiSequenceSHARKStarNodeChild;
                MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(confLowerBound, confUpperBound);
                child.setConfBounds(confBounds,
                        bound.sequence,
                        "drill-down");
                /*
                child.setPfuncBounds(confBounds.boltzmannWeight(msBound.bc),
                        bound.sequence,
                        "drill-down");
                 */
                child.setErrorBound(
                        msBound.bc.freeEnergy(confBounds.boltzmannWeight(msBound.bc).size(msBound.bc.mathContext)),
                        bound.sequence);
                //MultiSequenceSHARKStarNodeChild.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound,
                        //confUpperBound, msBound.bc.calc(confLowerBound), msBound.bc.calc(confUpperBound), bound.sequence, historyString);
                //MultiSequenceSHARKStarNodeChild.setBoundsFromConfLowerAndUpperWithHistory(confLowerBound, confUpperBound, msBound.bc.calc(confLowerBound), msBound.bc.calc(confUpperBound), bound.sequence, historyString);
                if (Double.isNaN(child.getPartialConfUpperBound()))
                    System.out.println("Huh!?");
                if (confLowerBound < bestChildLower) {
                    bestChild = child;
                }
                // collect the possible children
                if (child.getConfLowerBound(bound.sequence) < 0) {
                    children.add(child);
                }
                newNodes.add(child);

            }
            return bestChild;
        }
    }

    protected void boundLowestBoundConfUnderNode(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarBound msBound,
                                                 MultiSequenceSHARKStarNode startNode,
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
            ConfIndex<PartialConfAStarNode> index = new ConfIndex(RCs.getNumPos());
            curNode.index(index);

            if (curNode.getLevel() < RCs.getNumPos()) {
                MultiSequenceSHARKStarNode nextNode = drillDown(bound, msBound, newNodes, curNode);
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
        //addFullMinimizationsToCorrectionMatrix();
        return this.correctionMatrix;
    }

    /**
     * Takes the full minimizations from this tree, insert them into the correction matrix
     */
    private static void addFullMinimizationsToCorrectionMatrix(MultiSequenceSHARKStarBound msBound, SingleSequenceSHARKStarBound seqBound){
        for (MultiSequenceSHARKStarNode node : seqBound.finishedNodes){
            double nodeLB = msBound.minimizingEmat.confE(node.assignments);
            double nodeUB = msBound.rigidEmat.confE(node.assignments);
            double correctionSize = node.getPartialConfLowerBound() - nodeLB; //since the node is minimized, the first term is the minimized E
            if(correctionSize > 0 && nodeLB + correctionSize < nodeUB){
                msBound.correctionMatrix.setHigherOrder(new RCTuple(node.assignments), correctionSize);
            }
        }
        //captureSubtreeFullMinimizations(this.rootNode);
    }

    public List<TupE> getCorrections() {
        return new ArrayList<>(correctionMatrix.getAllCorrections());
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
        MultiSequenceSHARKStarNode resultNode = null;
        double upperBound = Double.NaN;
        double lowerBound = Double.NaN;
        String historyString = "Error!!";

        public boolean isValid() {
            return resultNode != null && !Double.isNaN(upperBound) && !Double.isNaN(lowerBound);
        }
    }

    //TODO: Simplify and clean up this method
    protected void processPartialConfNode(SingleSequenceSHARKStarBound bound, MultiSequenceSHARKStarBound msBound,
                                          List<MultiSequenceSHARKStarNode> newNodes,
                                          MultiSequenceSHARKStarNode curNode) {
        RCs RCs = bound.seqRCs;
        //debugPrint("Processing "+curNode.toSeqString(bound.sequence));
        // which pos to expand next?
        try (ObjectPool.Checkout<ScoreContext> checkout = msBound.contexts.autoCheckout()) {
            ScoreContext context = checkout.get();

            curNode.index(context.index);
            int nextPos = msBound.order.getNextPos(context.index, RCs);
            assert (!context.index.isDefined(nextPos));
            assert (context.index.isUndefined(nextPos));
            // score child nodes with tasks (possibly in parallel)

            for (int nextRc : bound.seqRCs.get(nextPos)) {
                double resultingUpper = Double.MAX_VALUE;
                double resultingLower = -Double.MAX_VALUE;
                String historyString = "";

                // index the parent, make the child
                curNode.index(context.index);
                MultiSequenceSHARKStarNode child;
                if(curNode.hasExistingChild(nextRc)){
                    child = curNode.getExistingChild(nextRc);
                }else{
                    child = curNode.assign(nextPos, nextRc);
                    // defer adding child to after we can determine the error of the child
                }

                // Set up corrections
                double confCorrection = Double.NEGATIVE_INFINITY;
                if(doCorrections)
                    confCorrection = msBound.correctionMatrix.confE(child.assignments);

                // score the child node differentially against the parent node
                // for now, don't calc differentially
                child.index(context.index);
                //double gScoreLB = context.partialConfLowerBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                //double gScoreUB = context.partialConfUpperBoundScorer.calcDifferential(context.index, RCs, nextPos, nextRc);
                double gScoreLB = context.partialConfLowerBoundScorer.calc(context.index, RCs);
                double gScoreUB = context.partialConfUpperBoundScorer.calc(context.index, RCs);
                double hScoreLB = context.lowerBoundScorer.calc(context.index, RCs);
                double hScoreUB = context.upperBoundScorer.calc(context.index, RCs);

                //Check to make sure our correction is reasonable
                if(gScoreLB < confCorrection && confCorrection < gScoreUB) {
                    //TODO: Warn on overcorrection?
                    child.setPartialConfLowerAndUpper(confCorrection, gScoreUB);
                    recordCorrection(gScoreLB+ hScoreLB, confCorrection - gScoreLB);
                }else {
                    child.setPartialConfLowerAndUpper(gScoreLB, gScoreUB);
                    if(doCorrections && confCorrection - gScoreLB > 1e-9)
                        System.err.printf("WARNING: tried to overcorrect %s: [%.3f, %.3f] by %.3f%n",
                                child.confToString(),
                                gScoreLB + hScoreLB,
                                gScoreUB + hScoreUB,
                                confCorrection - gScoreLB
                                );
                }

                double confLowerBound = child.getPartialConfLowerBound() + hScoreLB;
                double confUpperbound = child.getPartialConfUpperBound() + hScoreUB;
                //double lowerbound = minimizingEmat.confE(child.assignments);
                if (debug) {
                    historyString = String.format("%s: previous lower bound (none), g score %f, hscore %f, f score %f corrected score %f, from %s",
                            curNode.confToString(), curNode.getConfLowerBound(bound.sequence), gScoreLB, hScoreLB, gScoreLB + hScoreLB, confCorrection, getStackTrace());
                }
                //progress.reportInternalNode(child.level, child.getPartialConfLowerBound(), confLowerBound, queue.size(), children.size(), bound.getSequenceEpsilon());
                resultingLower = confLowerBound;
                resultingUpper= confUpperbound;
                // Make the MSSHARK node

                child.index(context.index);
                SimpleConfSpace.Position designPos = msBound.confSpace.positions.get(nextPos);
                int nextDesignIndex = msBound.order.getNextPos(context.index, bound.seqRCs);
                SimpleConfSpace.Position nextDesignPos = null;
                if (nextDesignIndex >= 0)
                    nextDesignPos = msBound.confSpace.positions.get(nextDesignIndex);
                /*
                MultiSequenceSHARKStarNode newChild = curNode.makeOrUpdateChild(child,
                        bound.sequence, resultingLower, resultingUpper, msBound.bc.calc(resultingLower), msBound.bc.calc(resultingUpper), designPos, nextDesignPos, bound.seqRCs);

                 */
                MathTools.DoubleBounds confBounds = new MathTools.DoubleBounds(resultingLower, resultingUpper);
                child.setConfBounds(confBounds,
                        bound.sequence,
                        "expansion");
                /*
                child.setPfuncBounds(confBounds.boltzmannWeight(msBound.bc),
                        bound.sequence,
                        "expansion");
                 */
                child.setErrorBound(
                        msBound.bc.freeEnergy(confBounds.boltzmannWeight(msBound.bc).size(msBound.bc.mathContext)),
                        bound.sequence);
                //newChild.setBoundsFromConfLowerAndUpperWithHistory(resultingLower,
                        //resultingUpper, msBound.bc.calc(resultingLower), msBound.bc.calc(resultingUpper), bound.sequence, historyString);
                //System.out.println("Created new child "+MultiSequenceSHARKStarNodeChild.toSeqString(bound.sequence));
                // collect the possible children
                if (child.getErrorBound(bound.sequence) < skipCutoff ||
                        bound.sequence.equals(this.precomputedSequence) ||
                        !skipAddingToFringe
                ) {
                    curNode.addChild(child, bound.sequence);
                }else{
                    debugPrint(String.format("Skipping child addition of node with no error"));
                }

                /*
                if (child.isMinimized(bound.sequence)) {
                    //newChild.computeEpsilonErrorBounds(bound.sequence);
                    //bound.addFinishedNode(newChild);
                    System.err.println(String.format("WARNING: perfectly bounded internal node: Seq %s %s, state %s, node %s: energy [%.3f, %.3f]",
                            bound.sequence,
                            bound.seqRCs,
                            this.cachePattern,
                            Arrays.toString(child.assignments),
                            child.getConfLowerBound(bound.sequence),
                            child.getConfUpperBound(bound.sequence)
                            ));
                }
                 */

                newNodes.add(child);

                //curNode.updateSubtreeBounds(bound.sequence);
                //printTree(bound.sequence, curNode);
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

    private void printMinimizationOutput(MultiSequenceSHARKStarNode node, double newConfLower, double oldgscore, SingleSequenceSHARKStarBound seqBound) {
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

    private boolean hasPrunedPair(ConfIndex<PartialConfAStarNode> confIndex, int nextPos, int nextRc) {

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
        public ConfIndex<PartialConfAStarNode> index;
        public AStarScorer<PartialConfAStarNode> partialConfLowerBoundScorer;
        public AStarScorer<PartialConfAStarNode> lowerBoundScorer;
        public AStarScorer<PartialConfAStarNode> upperBoundScorer;
        public AStarScorer<PartialConfAStarNode> partialConfUpperBoundScorer;
        public ConfEnergyCalculator ecalc;
        public BatchCorrectionMinimizer batcher;
    }

    public interface ScorerFactory<T extends PartialConfAStarNode> {
        AStarScorer<T> make(EnergyMatrix emat);
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

    public void printTimePerSequence(){
        System.out.println(String.format("For ID: %s, %d expanded, %d minimized, %d partials",
                this.cachePattern,
                this.state.numExpansions,
                this.state.numEnergiedConfs,
                this.state.numPartialMinimizations
                ));
        this.state.secondsPerSeq.forEach(
                (seq, time) -> System.out.println(String.format("Time spent computing for %s: %1.3e seconds",
                        seq, time)));
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

        Map<Sequence, Double> secondsPerSeq = new HashMap<>();
    }

    public static Pair<Double[], Double> generate1DRepresentation(SingleSequenceSHARKStarBound bound, int numBins, double cutoff){
        BoltzmannCalculator calc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        class occEntry{
            double occupancy;
            int numberConfs;

            occEntry(MultiSequenceSHARKStarNode node){
                if(bound.getUpperBound().compareTo(BigDecimal.ZERO) > 0)
                    this.occupancy = calc.calc(node.getConfLowerBound(bound.sequence)).divide(bound.getUpperBound(), PartitionFunction.decimalPrecision).doubleValue();
                else{
                    occupancy = 0;
                }
                this.numberConfs = MultiSequenceSHARKStarNode.computeNumConformations(node, bound.seqRCs).intValue();
            }

            occEntry(double occupancy, int numberConfs){
                this.occupancy = occupancy;
                this.numberConfs = numberConfs;
            }
        }

        List<occEntry> sortedEntries = Stream.of(bound.internalQueue, bound.leafQueue, bound.finishedNodes)
                .flatMap(Collection::stream)
                .parallel()
                .map(occEntry::new)
                .filter((e) -> e.occupancy / e.numberConfs > cutoff) // filter out elements with low occupancy
                .flatMap((e) -> Collections.nCopies(e.numberConfs, new occEntry(e.occupancy / e.numberConfs, 1)).stream()) // expand internal nodes
                .sorted(Comparator.comparingDouble((e) -> -1 * e.occupancy)) // sort by occupancy
                .collect(Collectors.toList());

        // generate a binned representation
        /*
        Double[] binned;
        if(sortedEntries.size() > 0) {
            int finalNumBins = Math.min(numBins, sortedEntries.size());
            int entriesPerBin = Math.max(sortedEntries.size() / finalNumBins, 1);
            binned = IntStream.range(0, finalNumBins)
                    .mapToObj((i) -> sortedEntries.subList(i * entriesPerBin, (i + 1) * entriesPerBin)
                            .parallelStream()
                            .map((e) -> e.occupancy)
                            .reduce(0.0, Double::sum))
                    .toArray(Double[]::new);
        }else{
            binned = new Double[] {};
        }
         */


        //or, just generate the first n
        Double[] firstn;
        if(sortedEntries.size() > 0) {
            int finalNumEntries = Math.min(numBins, sortedEntries.size());
            while(sortedEntries.size() < finalNumEntries){
                sortedEntries.add(new occEntry(0,1));
            }
                firstn = sortedEntries.subList(0, finalNumEntries).stream().map((n) -> n.occupancy).toArray(Double[]::new);
        }else{
            firstn = new Double[] {};
        }


        // compute the entropy
        Double entropy = -1 * BoltzmannCalculator.constRT * sortedEntries.parallelStream()
                .map((e) -> e.occupancy* Math.log(e.occupancy))
                .reduce(0.0, Double::sum);

        Pair<Double[], Double> binnedAndEntropy = new Pair<>(firstn, entropy);
        return binnedAndEntropy;

    }

}
