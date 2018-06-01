package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.ewakstar.EWAKStarLimitedSequenceTrie;
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
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.UpperBoundCalculator;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.function.Function;

/** Adapted from BBKStar.java
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 */

public class EWAKStarBBKStar {

    // *sigh* Java makes this stuff so verbose to do...
    // Kotlin would make this so much easier
    public static class Settings {

        public static class Builder {

            /**
             * The number of best (by K* score) sequences to evaluate before finishing
             */
            private EWAKStarLimitedSequenceTrie allowedSeqs;

            public Builder setAllowedSeqs(EWAKStarLimitedSequenceTrie seqs) {
                allowedSeqs = seqs;
                return this;
            }

            public EWAKStarBBKStar.Settings build() {
                return new EWAKStarBBKStar.Settings(allowedSeqs);
            }
        }

        public final EWAKStarLimitedSequenceTrie allowedSeqs;

        public Settings(EWAKStarLimitedSequenceTrie seqs) {
            this.allowedSeqs = seqs;
        }
    }

    public class ConfSpaceInfo {

        public final KStar.ConfSpaceType type;
        public final SimpleConfSpace confSpace;
        public final ConfEnergyCalculator minimizingConfEcalc;
        public final ConfEnergyCalculator rigidConfEcalc;

        public EnergyMatrix rigidEmat = null;
        public final EnergyMatrix emat;

        public BigDecimal stabilityThreshold = null;

        public ConfSpaceInfo(KStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator minimizingConfEcalc, ConfEnergyCalculator rigidConfECalc, EnergyMatrix emat) {
            this.type = type;
            this.confSpace = confSpace;
            this.minimizingConfEcalc = minimizingConfEcalc;
            this.emat = emat;
            this.rigidConfEcalc = rigidConfECalc;
        }

        public void calcEmatsIfNeeded() {
            if (rigidEmat == null) {
                SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc);
                if (kstarSettings.energyMatrixCachePattern != null) {
                    builder.setCacheFile(new File(kstarSettings.applyEnergyMatrixCachePattern(type.name().toLowerCase() + ".rigidNegated")));
                }
                rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
                        .build()
                        .calcEnergyMatrix();

                // negate the rigid energy matrix, so we can convince A* to return scores
                // in descending order instead of ascending order
                // (of course, we'll have to negate the scores returned by A* too)
                //rigidEmat.negate();
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
        public int numCurSeqs = 0;

        /** for comparing in the tree, higher is first */
        public double score;

        /** signals whether or not partition function values are allowed the stability threshold */
        public boolean isUnboundUnstable;

        protected Node(Sequence sequence, EWAKStarBBKStar.ConfDBs confdbs) {

            if(sequence.isFullyAssigned()){
                numCurSeqs++;
            }
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

            // get the possible assignments// get the possible assignments
            //only allow sub-sequences that exist in our limited sequence space
            Set<String> resTypes;


            String subSeq = sequence.getAssignedResTypes();

            if(assignPos.index == 0 || assignPos.resFlex.resTypes.size() == 1)
                resTypes = new HashSet<>(assignPos.resFlex.resTypes);
            else {
                resTypes = filterOnPreviousSeqs();
            }
            // for each assignment...
            for (String resType : resTypes) {

                // update the sequence with this assignment
                Sequence s = sequence.makeMutatedSequence(assignPos, resType);

                if (s.isFullyAssigned()) {
                    // fully assigned, make single sequence node

                    children.add(new SingleSequenceNode(s, confdbs));

                } else if (kstarSettings.useExact && s.countMutations() == kstarSettings.maxSimultaneousMutations) {
                    // mutation limit reached, fill unassigned positions with wild-type
                    s.fillWildType();
                    children.add(new EWAKStarBBKStar.SingleSequenceNode(s, confdbs));
                } else {
                    // still partial sequence, make multi-sequence node
                    children.add(new EWAKStarBBKStar.MultiSequenceNode(s, confdbs));
                }
            }

            return children;
        }

        private Set<String> filterOnPreviousSeqs(){

            String subSeq = sequence.getAssignedResTypes();
            //"PHE PHE ILE THR PHE ASP GLU THR"
            //WT: "PHE LYS ILE THR PHE ASP GLU THR"
            Set<String> resTypes = bbkstarSettings.allowedSeqs.getSeq(subSeq);
            return resTypes;
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
                    EWAKStarPartitionFunction.decimalPrecision
            ).doubleValue();
            isUnboundUnstable = false;
        }

        private BigDecimal calcLowerBound(EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute lower bounds on pfuncs, we'll do the usual lower bound calculation,
            // but use rigid energies instead of minimized energies

            RCs rcs = sequence.makeRCs();
            ConfSearch astar = new ConfAStarTree.Builder(info.rigidEmat, rcs)
                    .setTraditional()
                    .build();
            BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

            BigMath m = new BigMath(PartitionFunction.decimalPrecision)
                    .set(0.0);
            for (int i=0; i<numConfs; i++) {
                ConfSearch.ScoredConf conf = astar.nextConf();
                if (conf == null) {
                    break;
                }

                m.add(bcalc.calc(conf.getScore()));
            }
            return m.get();
        }

        private BigDecimal calcUpperBound(EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute upper bounds on pfuncs,
            // we'll use the upper bound calculator in the usual way

            RCs rcs = sequence.makeRCs();
            ConfSearch astar = new ConfAStarTree.Builder(info.emat, rcs)
                    .setTraditional()
                    .build();
            UpperBoundCalculator calc = new UpperBoundCalculator(
                    astar,
                    rcs.getNumConformations()
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

        public final EWAKStarPartitionFunction protein;
        public final EWAKStarPartitionFunction ligand;
        public final EWAKStarPartitionFunction complex;

        public SingleSequenceNode(Sequence sequence, EWAKStarBBKStar.ConfDBs confdbs) {
            super(sequence, confdbs);

            // make the partition functions
            this.protein = makePfunc(proteinPfuncs, EWAKStarBBKStar.this.protein, proteinSequence, confdbs.protein, ematP);
            this.ligand = makePfunc(ligandPfuncs, EWAKStarBBKStar.this.ligand, ligandSequence, confdbs.ligand, ematL);
            this.complex = makePfunc(complexPfuncs, EWAKStarBBKStar.this.complex, sequence, confdbs.complex, ematPL);
        }

        private EWAKStarPartitionFunction makePfunc(Map<Sequence,EWAKStarPartitionFunction> pfuncCache, EWAKStarBBKStar.ConfSpaceInfo info, Sequence sequence, ConfDB confdb, EnergyMatrix curEmat) {

            // first check the cache
            EWAKStarPartitionFunction pfunc = pfuncCache.get(sequence);
            if (pfunc != null) {
                return pfunc;
            }

            // cache miss, need to compute the partition function

            // make the partition function
            ConfSearch astar = confSearchFactory.make(curEmat, sequence.makeRCs());
            EWAKStarGradientDescentPfunc newPfunc = new EWAKStarGradientDescentPfunc(astar, info.minimizingConfEcalc);

            newPfunc.setReportProgress(kstarSettings.showPfuncProgress);
            newPfunc.init(kstarSettings.eW, kstarSettings.epsilon, astar.getNumConformations());
            pfuncCache.put(sequence, newPfunc);
            return newPfunc;
        }

        @Override
        public void estimateScore() {

            // tank the sequence if either unbound strand is unstable
            // yeah, we haven't refined any pfuncs yet this estimation,
            // but since pfuncs get cached, check before we do any more estimation
            if (protein.getStatus() == EWAKStarPartitionFunction.Status.Unstable
                    || ligand.getStatus() == EWAKStarPartitionFunction.Status.Unstable) {
                score = Double.NEGATIVE_INFINITY;
                isUnboundUnstable = true;
                return;
            }

            // refine the pfuncs if needed
            if (protein.getStatus().canContinue()) {
                protein.compute(kstarSettings.maxPFConfs);

                // tank the sequence if the unbound protein is unstable
                if (protein.getStatus() == EWAKStarPartitionFunction.Status.Unstable) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            if (ligand.getStatus().canContinue()) {
                ligand.compute(kstarSettings.maxPFConfs);

                // tank the sequence if the unbound ligand is unstable
                if (ligand.getStatus() == EWAKStarPartitionFunction.Status.Unstable) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            if (complex.getStatus().canContinue()) {
                complex.compute(kstarSettings.maxPFConfs);
            }

            // update the score
            score = Math.log10(makeKStarScore().upperBound.doubleValue());
            isUnboundUnstable = false;
        }

        public EWAKStarScore computeScore() {

            // refine the pfuncs until done
            while (protein.getStatus().canContinue()) {
                protein.compute(kstarSettings.maxPFConfs);
            }
            while (ligand.getStatus().canContinue()) {
                ligand.compute(kstarSettings.maxPFConfs);
            }
            while (complex.getStatus().canContinue()) {
                complex.compute(kstarSettings.maxPFConfs);
            }

            // update the score
            EWAKStarScore kstarScore = makeKStarScore();
            score = Math.log10(kstarScore.upperBound.doubleValue());
            return kstarScore;
        }

        public EWAKStarScore makeKStarScore() {
            return new EWAKStarScore(protein.makeResult(), ligand.makeResult(), complex.makeResult());
        }

        public EWAKStarBBKStar.PfuncsStatus getStatus() {

            // aggregate pfunc statuses
            if (protein.getStatus() == EWAKStarPartitionFunction.Status.Estimated
                    && ligand.getStatus() == EWAKStarPartitionFunction.Status.Estimated
                    && complex.getStatus() == EWAKStarPartitionFunction.Status.Estimated) {
                return EWAKStarBBKStar.PfuncsStatus.Estimated;
            } else if (protein.getStatus() == EWAKStarPartitionFunction.Status.Estimating
                    || ligand.getStatus() == EWAKStarPartitionFunction.Status.Estimating
                    || complex.getStatus() == EWAKStarPartitionFunction.Status.Estimating) {
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

    /** Calculates the continuous-rotamer (minimized) energy for a molecule */
    public final EnergyCalculator minimizingEcalc;

    /** Calculates the continuous-rotamer (minimized) energy for a molecule */
    public final EnergyMatrix ematPL;

    /** Calculates the continuous-rotamer (minimized) energy for a molecule */
    public final EnergyMatrix ematP;

    /** Calculates the continuous-rotamer (minimized) energy for a molecule */
    public final EnergyMatrix ematL;

    /** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
    public final EWAKStar.ConfSearchFactory confSearchFactory;

    /** Optional and overridable settings for BBK*, shared with K* */
    public final EWAKStar.Settings kstarSettings;
    public final EWAKStarBBKStar.Settings bbkstarSettings;

    public int maxNumSeqs = 0;

    // TODO: caching these will keep lots of A* trees in memory. is that a problem?
    private final Map<Sequence,EWAKStarPartitionFunction> proteinPfuncs;
    private final Map<Sequence,EWAKStarPartitionFunction> ligandPfuncs;
    private final Map<Sequence,EWAKStarPartitionFunction> complexPfuncs;

	public EWAKStarBBKStar(int maxNumSeqs, SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator minimizingEcalc, ConfEnergyCalculator confECalcPL, ConfEnergyCalculator confECalcP, ConfEnergyCalculator confECalcL, ConfEnergyCalculator confRigidECalcPL, ConfEnergyCalculator confRigidECalcP, ConfEnergyCalculator confRigidECalcL,  EWAKStar.ConfSearchFactory confSearchFactory, EWAKStarBBKStar.Settings bbkstarSettings, EWAKStar.Settings kstarSettings, EnergyMatrix ematPL, EnergyMatrix ematP, EnergyMatrix ematL) {

        this.protein = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Protein,
                protein,
                confECalcP,
                confRigidECalcP,
                ematP
        );
        this.ligand = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Ligand,
                ligand,
                confECalcL,
                confRigidECalcL,
                ematL
        );
        this.complex = new EWAKStarBBKStar.ConfSpaceInfo(
                KStar.ConfSpaceType.Complex,
                complex,
                confECalcPL,
                confRigidECalcPL,
                ematPL
        );

        this.maxNumSeqs = maxNumSeqs;
        this.bbkstarSettings = bbkstarSettings;
        this.minimizingEcalc = minimizingEcalc;
        this.confSearchFactory = confSearchFactory;
        this.kstarSettings = kstarSettings;
        this.ematPL = ematPL;
        this.ematL = ematL;
        this.ematP = ematP;

        proteinPfuncs = new HashMap<>();
        ligandPfuncs = new HashMap<>();
        complexPfuncs = new HashMap<>();
    }

    public List<EWAKStar.ScoredSequence> run() {

        protein.calcEmatsIfNeeded();
        ligand.calcEmatsIfNeeded();
        complex.calcEmatsIfNeeded();

        List<EWAKStar.ScoredSequence> scoredSequences = new ArrayList<>();

        ConfDBs confdbs = new EWAKStarBBKStar.ConfDBs();

        // start the BBK* tree with the root node
        PriorityQueue<EWAKStarBBKStar.Node> tree = new PriorityQueue<>();
        tree.add(new EWAKStarBBKStar.MultiSequenceNode(complex.confSpace.makeUnassignedSequence(), confdbs));

        // start searching the tree
        System.out.println("Computing K* scores for the best sequences to within an energy window of " + kstarSettings.eW + " kcal, with a max of "+kstarSettings.maxPFConfs+ " conformations, and an epsilon of "+kstarSettings.epsilon+"...");
        kstarSettings.scoreWriters.writeHeader();
        Boolean wtSeqFound = false;
        if(kstarSettings.wtBenchmark) {
            while (!tree.isEmpty() && !wtSeqFound) {

                // get the next node
                EWAKStarBBKStar.Node node = tree.poll();

                if (node.isWTSeq) {
                    wtSeqFound = true;
                }

                if (node instanceof EWAKStarBBKStar.SingleSequenceNode) {
                    EWAKStarBBKStar.SingleSequenceNode ssnode = (EWAKStarBBKStar.SingleSequenceNode) node;

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
                    EWAKStarBBKStar.MultiSequenceNode msnode = (EWAKStarBBKStar.MultiSequenceNode) node;

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
        } else {
            while (!tree.isEmpty() && scoredSequences.size() < maxNumSeqs) {

                // get the next node
                EWAKStarBBKStar.Node node = tree.poll();

                if (node instanceof EWAKStarBBKStar.SingleSequenceNode) {
                    EWAKStarBBKStar.SingleSequenceNode ssnode = (EWAKStarBBKStar.SingleSequenceNode) node;

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
                    EWAKStarBBKStar.MultiSequenceNode msnode = (EWAKStarBBKStar.MultiSequenceNode) node;

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
        }


        if (tree.isEmpty()) {
            // all is well, we just don't have that many sequences in the design
            System.out.println("All " + scoredSequences.size() + " sequences calculated. EWAK* complete.");
        } else if (wtSeqFound && scoredSequences.size()==1) {
            System.out.println("No K* scores found that are better than the wild-type sequence.");
        } else if (wtSeqFound) {
                System.out.println("Found K* score estimates for all " + scoredSequences.size() + " sequences with scores greater than that of the wild-type sequence.");
        } else if (scoredSequences.size() >= maxNumSeqs) {
            System.out.println("Found K* score estimates for top "+ maxNumSeqs +" sequences.");
        } else
            throw new Error("EWAK* ended, but the tree isn't empty and we didn't return all of the sequences. This is a bug.");

        return scoredSequences;
    }

    private void reportSequence(EWAKStarBBKStar.SingleSequenceNode ssnode, List<EWAKStar.ScoredSequence> scoredSequences) {

        EWAKStarScore kstarScore = ssnode.makeKStarScore();
        scoredSequences.add(new EWAKStar.ScoredSequence(ssnode.sequence, kstarScore));

        kstarSettings.scoreWriters.writeScore(new EWAKStarScoreWriter.ScoreInfo(
                scoredSequences.size() - 1,
                0,
                ssnode.sequence,
                complex.confSpace,
                kstarScore
        ));
    }
}