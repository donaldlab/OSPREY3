package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.UpperBoundCalculator;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.function.Function;

/** Adapted from BBKStar.java by lowegard **/

public class EWAKStarBBKStar {

    // *sigh* Java makes this stuff so verbose to do...
    // Kotlin would make this so much easier
    public static class Settings {

        public static class Builder {

            /** The number of best (by K* score) sequences to evaluate before finishing */
            private double ewForKStarScoreCalc = 2.0;

            /**
             * The number of conformations to evaluate per batch of partition function refinement
             *
             * For best results, make this a multiple of the available parallelism. e.g., if using 4 threads,
             * try a batch size of 4, 8, 12, 16, etc.
             */
            private int numConfsPerBatch = 8;

            public EWAKStarBBKStar.Settings.Builder setNumConfsPerBatch(int val) {
                numConfsPerBatch = val;
                return this;
            }

            public EWAKStarBBKStar.Settings build() {
                return new EWAKStarBBKStar.Settings(ewForKStarScoreCalc, numConfsPerBatch);
            }
        }

        public final double ewForKStarScoreCalc;
        public final int numConfsPerBatch;

        public Settings(double ewForKStarScoreCalc, int numConfsPerBatch) {
            this.ewForKStarScoreCalc = ewForKStarScoreCalc;
            this.numConfsPerBatch = numConfsPerBatch;
        }
    }

    public class ConfSpaceInfo {

        public final KStar.ConfSpaceType type;
        public final SimpleConfSpace confSpace;
        public final ConfEnergyCalculator rigidConfEcalc;
        public final ConfEnergyCalculator minimizingConfEcalc;

        public EnergyMatrix rigidNegatedEmat = null;
        public EnergyMatrix minimizedEmat = null;
        public BigDecimal stabilityThreshold = null;

        public ConfSpaceInfo(KStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator rigidConfEcalc, ConfEnergyCalculator minimizingConfEcalc) {
            this.type = type;
            this.confSpace = confSpace;
            this.rigidConfEcalc = rigidConfEcalc;
            this.minimizingConfEcalc = minimizingConfEcalc;
        }

        public void calcEmatsIfNeeded() {
            if (rigidNegatedEmat == null) {
                SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc);
                if (kstarSettings.energyMatrixCachePattern != null) {
                    builder.setCacheFile(new File(kstarSettings.applyEnergyMatrixCachePattern(type.name().toLowerCase() + ".rigidNegated")));
                }
                rigidNegatedEmat = builder.build().calcEnergyMatrix();

                // negate the rigid energy matrix, so we can convince A* to return scores
                // in descending order instead of ascending order
                // (of course, we'll have to negate the scores returned by A* too)
                rigidNegatedEmat.negate();
            }
            if (minimizedEmat == null) {
                SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc);
                if (kstarSettings.energyMatrixCachePattern != null) {
                    builder.setCacheFile(new File(kstarSettings.applyEnergyMatrixCachePattern(type.name().toLowerCase())));
                }
                minimizedEmat = builder.build().calcEnergyMatrix();
            }
        }
    }

    private class ConfDBs {
        public ConfDB protein = null;
        public ConfDB ligand = null;
        public ConfDB complex = null;
    }

    public static enum PfuncsStatus {
        Estimating,
        Estimated,
        Blocked
    }

    private abstract class Node implements Comparable<EWAKStarBBKStar.Node> {

        public final Sequence sequence;
        public final Sequence proteinSequence;
        public final Sequence ligandSequence;
        public final EWAKStarBBKStar.ConfDBs confdbs;
        public final boolean isWTSeq;

        /** for comparing in the tree, higher is first */
        public double score;

        /** signals whether or not partition function values are allowed the stability threshold */
        public boolean isUnboundUnstable;

        protected Node(Sequence sequence, EWAKStarBBKStar.ConfDBs confdbs) {

            this.sequence = sequence;
            this.isWTSeq = sequence.toString().equals(EWAKStarBBKStar.this.complex.confSpace.makeWildTypeSequence().toString());
            this.confdbs = confdbs;

            // split complex sequence into protein/ligand sequences
            proteinSequence = sequence.filter(EWAKStarBBKStar.this.protein.confSpace);
            ligandSequence = sequence.filter(EWAKStarBBKStar.this.ligand.confSpace);

            this.score = 0;
        }

        @Override
        public int compareTo(EWAKStarBBKStar.Node other) {
            // negate for descending sort
            return -Double.compare(this.score, other.score);
        }

        public abstract void estimateScore();
    }

    public class MultiSequenceNode extends EWAKStarBBKStar.Node {

        public MultiSequenceNode(Sequence sequence, EWAKStarBBKStar.ConfDBs confdbs) {
            super(sequence, confdbs);
        }

        public List<EWAKStarBBKStar.Node> makeChildren() {

            List<EWAKStarBBKStar.Node> children = new ArrayList<>();

            // pick the next design position
            // TODO: dynamic A*?
            List<SimpleConfSpace.Position> positions = complex.confSpace.positions;
            SimpleConfSpace.Position assignPos = positions.stream()
                    .filter((pos) -> !sequence.isAssigned(pos))
                    .findFirst()
                    .orElseThrow(() -> new IllegalStateException("no design positions left to choose"));

            // get the possible assignments
            Set<String> resTypes = new HashSet<>(assignPos.resFlex.resTypes);

            // for each assignment...
            for (String resType : resTypes) {

                // update the sequence with this assignment
                Sequence s = sequence.makeMutatedSequence(assignPos, resType);

                if (s.isFullyAssigned()) {

                    // fully assigned, make single sequence node
                    children.add(new SingleSequenceNode(s, confdbs));

                } else {

                    // still partial sequence, make multi-sequence node
                    children.add(new EWAKStarBBKStar.MultiSequenceNode(s, confdbs));
                }
            }

            return children;
        }

        @Override
        public void estimateScore() {

            // TODO: expose setting?
            // NOTE: for the correctness of the bounds, the number of confs must be the same for every node
            // meaning, it might not be sound to do epsilon-based iterative approximations here
            final int numConfs = 1000;

            if (protein.stabilityThreshold != null) {

                // tank the sequence if the protein is unstable
                BigDecimal proteinUpperBound = calcUpperBound(protein, proteinSequence, numConfs);
                if (MathTools.isLessThan(proteinUpperBound, protein.stabilityThreshold)) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            BigDecimal proteinLowerBound = calcLowerBound(protein, proteinSequence, numConfs);

            // if the first few conf upper bound scores (for the pfunc lower bound) are too high,
            // then the K* upper bound is also too high
            if (MathTools.isZero(proteinLowerBound)) {
                score = Double.POSITIVE_INFINITY;
                isUnboundUnstable = false;
                return;
            }

            if (ligand.stabilityThreshold != null) {

                // tank the sequence if the ligand is unstable
                BigDecimal ligandUpperBound = calcUpperBound(ligand, ligandSequence, numConfs);
                if (MathTools.isLessThan(ligandUpperBound, ligand.stabilityThreshold)) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            BigDecimal ligandLowerBound = calcLowerBound(ligand, ligandSequence, numConfs);

            // if the first few conf upper bound scores (for the pfunc lower bound) are too high,
            // then the K* upper bound is also too high
            if (MathTools.isZero(ligandLowerBound)) {
                score = Double.POSITIVE_INFINITY;
                isUnboundUnstable = false;
                return;
            }

            BigDecimal complexUpperBound = calcUpperBound(complex, sequence, numConfs);

            // compute the node score
            score = MathTools.bigDivideDivide(
                    complexUpperBound,
                    proteinLowerBound,
                    ligandLowerBound,
                    PartitionFunction.decimalPrecision
            ).doubleValue();
            isUnboundUnstable = false;
        }

        private BigDecimal calcLowerBound(EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute lower bounds on pfuncs, we'll use the upper bound calculator,
            // but feed it upper-bound scores in order of greatest upper bound
            // instead of the usual lower-bound scores in order of smallest lower bound
            // which means we need to feed A* with a rigid, negated energy matrix
            // and then negate all the scores returned by A*
            Function<ConfSearch,ConfSearch> astarNegater = (confSearch) -> new ConfSearch() {
                @Override
                public ScoredConf nextConf() {
                    ScoredConf conf = confSearch.nextConf();
                    conf.setScore(-conf.getScore());
                    return conf;
                }

                @Override
                public BigInteger getNumConformations() {
                    return confSearch.getNumConformations();
                }

                @Override
                public List<ScoredConf> nextConfs(double maxEnergy) {
                    throw new UnsupportedOperationException("some lazy programmer didn't implement this =P");
                }
            };

            UpperBoundCalculator calc = new UpperBoundCalculator(
                    astarNegater.apply(confSearchFactory.make(info.rigidNegatedEmat, sequence.makeRCs()))
            );
            calc.run(numConfs);
            return calc.totalBound;
        }

        private BigDecimal calcUpperBound(EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute upper bounds on pfuncs,
            // we'll use the upper bound calculator in the usual way

            UpperBoundCalculator calc = new UpperBoundCalculator(
                    confSearchFactory.make(info.minimizedEmat, sequence.makeRCs())
            );
            calc.run(numConfs);
            return calc.totalBound;
        }

        @Override
        public String toString() {
            return String.format("MultiSequenceNode[score=%12.6f, seq=%s]",
                    score,
                    sequence
            );
        }
    }

    public class SingleSequenceNode extends EWAKStarBBKStar.Node {

        public final PartitionFunction protein;
        public final PartitionFunction ligand;
        public final PartitionFunction complex;

        public SingleSequenceNode(Sequence sequence, EWAKStarBBKStar.ConfDBs confdbs) {
            super(sequence, confdbs);

            // make the partition functions
            this.protein = makePfunc(proteinPfuncs, EWAKStarBBKStar.this.protein, proteinSequence, confdbs.protein);
            this.ligand = makePfunc(ligandPfuncs, EWAKStarBBKStar.this.ligand, ligandSequence, confdbs.ligand);
            this.complex = makePfunc(complexPfuncs, EWAKStarBBKStar.this.complex, sequence, confdbs.complex);
        }

        private PartitionFunction makePfunc(Map<Sequence,PartitionFunction> pfuncCache, EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, ConfDB confdb) {

            // first check the cache
            PartitionFunction pfunc = pfuncCache.get(sequence);
            if (pfunc != null) {
                return pfunc;
            }

            // cache miss, need to compute the partition function

            // make the partition function
            GradientDescentPfunc gdpfunc = new GradientDescentPfunc(confSearchFactory.make(info.minimizedEmat, sequence.makeRCs()), info.minimizingConfEcalc);
            gdpfunc.setReportProgress(kstarSettings.showPfuncProgress);
            if (confdb != null) {
                gdpfunc.setConfTable(confdb.getSequence(sequence));
            }
            gdpfunc.init(kstarSettings.epsilon, info.stabilityThreshold);
            pfuncCache.put(sequence, gdpfunc);
            return gdpfunc;
        }

        @Override
        public void estimateScore() {

            // tank the sequence if either unbound strand is unstable
            // yeah, we haven't refined any pfuncs yet this estimation,
            // but since pfuncs get cached, check before we do any more estimation
            if (protein.getStatus() == PartitionFunction.Status.Unstable
                    || ligand.getStatus() == PartitionFunction.Status.Unstable) {
                score = Double.NEGATIVE_INFINITY;
                isUnboundUnstable = true;
                return;
            }

            // refine the pfuncs if needed
            if (protein.getStatus().canContinue()) {
                protein.compute(bbkstarSettings.numConfsPerBatch);

                // tank the sequence if the unbound protein is unstable
                if (protein.getStatus() == PartitionFunction.Status.Unstable) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            if (ligand.getStatus().canContinue()) {
                ligand.compute(bbkstarSettings.numConfsPerBatch);

                // tank the sequence if the unbound ligand is unstable
                if (ligand.getStatus() == PartitionFunction.Status.Unstable) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            if (complex.getStatus().canContinue()) {
                complex.compute(bbkstarSettings.numConfsPerBatch);
            }

            // update the score
            score = Math.log10(makeKStarScore().upperBound.doubleValue());
            isUnboundUnstable = false;
        }

        public KStarScore computeScore() {

            // refine the pfuncs until done
            while (protein.getStatus().canContinue()) {
                protein.compute(bbkstarSettings.numConfsPerBatch);
            }
            while (ligand.getStatus().canContinue()) {
                ligand.compute(bbkstarSettings.numConfsPerBatch);
            }
            while (complex.getStatus().canContinue()) {
                complex.compute(bbkstarSettings.numConfsPerBatch);
            }

            // update the score
            KStarScore kstarScore = makeKStarScore();
            score = Math.log10(kstarScore.upperBound.doubleValue());
            return kstarScore;
        }

        public KStarScore makeKStarScore() {
            return new KStarScore(protein.makeResult(), ligand.makeResult(), complex.makeResult());
        }

        public EWAKStarBBKStar.PfuncsStatus getStatus() {

            // aggregate pfunc statuses
            if (protein.getStatus() == PartitionFunction.Status.Estimated
                    && ligand.getStatus() == PartitionFunction.Status.Estimated
                    && complex.getStatus() == PartitionFunction.Status.Estimated) {
                return EWAKStarBBKStar.PfuncsStatus.Estimated;
            } else if (protein.getStatus() == PartitionFunction.Status.Estimating
                    || ligand.getStatus() == PartitionFunction.Status.Estimating
                    || complex.getStatus() == PartitionFunction.Status.Estimating) {
                return EWAKStarBBKStar.PfuncsStatus.Estimating;
            } else {
                return EWAKStarBBKStar.PfuncsStatus.Blocked;
            }
        }

        @Override
        public String toString() {
            return String.format("SingleSequenceNode[score=%12.6f, seq=%s, K*=%s]",
                    score,
                    sequence,
                    makeKStarScore()
            );
        }
    }

    /** A configuration space containing just the protein strand */
    public final EWAKStarBBKStar.ConfSpaceInfo protein;

    /** A configuration space containing just the ligand strand */
    public final EWAKStarBBKStar.ConfSpaceInfo ligand;

    /** A configuration space containing both the protein and ligand strands */
    public final EWAKStarBBKStar.ConfSpaceInfo complex;

    /** Calculates the rigid-rotamer energy for a molecule */
    public final EnergyCalculator rigidEcalc;

    /** Calculates the continuous-rotamer (minimized) energy for a molecule */
    public final EnergyCalculator minimizingEcalc;

    /** A function that makes a ConfEnergyCalculator with the desired options */
    public final EWAKStar.ConfEnergyCalculatorFactory confEcalcFactory;

    /** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
    public final EWAKStar.ConfSearchFactory confSearchFactory;

    /** Optional and overridable settings for BBK*, shared with K* */
    public final EWAKStar.Settings kstarSettings;

    /** Optional and overridable settings for BBK* */
    public final EWAKStarBBKStar.Settings bbkstarSettings;

    // TODO: caching these will keep lots of A* trees in memory. is that a problem?
    private final Map<Sequence,PartitionFunction> proteinPfuncs;
    private final Map<Sequence,PartitionFunction> ligandPfuncs;
    private final Map<Sequence,PartitionFunction> complexPfuncs;

	public EWAKStarBBKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator rigidEcalc, EnergyCalculator minimizingEcalc, EWAKStar.ConfEnergyCalculatorFactory confEcalcFactory, EWAKStar.ConfSearchFactory confSearchFactory, EWAKStar.Settings kstarSettings, EWAKStarBBKStar.Settings bbkstarSettings) {

        this.protein = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Protein,
                protein,
                confEcalcFactory.make(protein, rigidEcalc),
                confEcalcFactory.make(protein, minimizingEcalc)
        );
        this.ligand = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Ligand,
                ligand,
                confEcalcFactory.make(ligand, rigidEcalc),
                confEcalcFactory.make(ligand, minimizingEcalc)
        );
        this.complex = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Complex,
                complex,
                confEcalcFactory.make(complex, rigidEcalc),
                confEcalcFactory.make(complex, minimizingEcalc)
        );
        this.rigidEcalc = rigidEcalc;
        this.minimizingEcalc = minimizingEcalc;
        this.confEcalcFactory = confEcalcFactory;
        this.confSearchFactory = confSearchFactory;
        this.kstarSettings = kstarSettings;
        this.bbkstarSettings = bbkstarSettings;

        proteinPfuncs = new HashMap<>();
        ligandPfuncs = new HashMap<>();
        complexPfuncs = new HashMap<>();
    }

    public List<EWAKStar.ScoredSequence> run(ArrayList<Sequence> inputSeqs) {

        protein.calcEmatsIfNeeded();
        ligand.calcEmatsIfNeeded();
        complex.calcEmatsIfNeeded();

        List<EWAKStar.ScoredSequence> scoredSequences = new ArrayList<>();

        ConfDBs confdbs = new EWAKStarBBKStar.ConfDBs();

        // calculate wild-type first
        System.out.println("computing K* score for the wild-type sequence...");
        EWAKStarBBKStar.SingleSequenceNode wildTypeNode = new EWAKStarBBKStar.SingleSequenceNode(Sequence.makeWildType(complex.confSpace), confdbs);
        KStarScore wildTypeScore = wildTypeNode.computeScore();
        kstarSettings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
                -1,
                0,
                wildTypeNode.sequence,
                complex.confSpace,
                wildTypeScore
        ));

        if (kstarSettings.stabilityThreshold != null) {
            BigDecimal stabilityThresholdFactor = new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(kstarSettings.stabilityThreshold);
            protein.stabilityThreshold = wildTypeScore.protein.values.calcLowerBound().multiply(stabilityThresholdFactor);
            ligand.stabilityThreshold = wildTypeScore.ligand.values.calcLowerBound().multiply(stabilityThresholdFactor);
        }

        // start the BBK* tree with the root node
        PriorityQueue<EWAKStarBBKStar.Node> tree = new PriorityQueue<>();
        tree.add(new EWAKStarBBKStar.MultiSequenceNode(complex.confSpace.makeUnassignedSequence(), confdbs));

        // start searching the tree
        System.out.println("computing K* scores for the best sequences to within an energy window of " + kstarSettings.eW + "  kcal...");
        kstarSettings.scoreWriters.writeHeader();
        Boolean wtSeqFound = false;
        while (!tree.isEmpty() && !wtSeqFound) {

            // get the next node
            EWAKStarBBKStar.Node node = tree.poll();

            if (node.isWTSeq){
                wtSeqFound = true;
            }

            if (node instanceof EWAKStarBBKStar.SingleSequenceNode) {
                EWAKStarBBKStar.SingleSequenceNode ssnode = (EWAKStarBBKStar.SingleSequenceNode)node;

                // single-sequence node
                switch (ssnode.getStatus()) {
                    case Estimated:

                        // sequence is finished, return it!
                        reportSequence(ssnode, scoredSequences);

                        break;
                    case Estimating:

                        // needs more estimation, catch-and-release
                        ssnode.estimateScore();
                        if (!ssnode.isUnboundUnstable) {
                            tree.add(ssnode);
                        }

                        break;
                    case Blocked:

                        // from here on out, it's all blocked sequences
                        // so it's ok to put them in the sorted order now
                        reportSequence(ssnode, scoredSequences);
                }

            } else if (node instanceof EWAKStarBBKStar.MultiSequenceNode) {
                EWAKStarBBKStar.MultiSequenceNode msnode = (EWAKStarBBKStar.MultiSequenceNode)node;

                // partial sequence, expand children
                // TODO: parallelize the multi-sequence node scoring here?
                for (EWAKStarBBKStar.Node child : msnode.makeChildren()) {
                    child.estimateScore();
                    if (!child.isUnboundUnstable) {
                        tree.add(child);
                    }
                }
            }
        }


        if (tree.isEmpty()) {
            // all is well, we just don't have that many sequences in the design
            System.out.println("Tried to find all sequences, but design flexibility and sequence filters only allowed " + scoredSequences.size() + " sequences.");
        } else {
            throw new Error("BBK* ended, but the tree isn't empty and we didn't return enough sequences. This is a bug.");
        }

        return scoredSequences;
    }

    private void reportSequence(EWAKStarBBKStar.SingleSequenceNode ssnode, List<EWAKStar.ScoredSequence> scoredSequences) {

        KStarScore kstarScore = ssnode.makeKStarScore();
        scoredSequences.add(new EWAKStar.ScoredSequence(ssnode.sequence, kstarScore));

        kstarSettings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
                scoredSequences.size() - 1,
                0,
                ssnode.sequence,
                complex.confSpace,
                kstarScore
        ));
    }
}