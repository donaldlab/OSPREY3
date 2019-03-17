/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfSearchCache;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.UpperBoundCalculator;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;
import java.util.function.Function;

/** Adapted from BBKStar.java by lowegard **/

public class EWAKStarBBKStar {

    // *sigh* Java makes this stuff so verbose to do...
    // Kotlin would make this so much easier
    public static class Settings {

        public static class Builder {

            /** The number of best (by K* score) sequences to evaluate before finishing */
            private int numBestSequences = 1;
            private EwakstarLimitedSequenceTrie allowedSeqs;

            /**
             * The number of conformations to evaluate per batch of partition function refinement
             *
             * For best results, make this a multiple of the available parallelism. e.g., if using 4 threads,
             * try a batch size of 4, 8, 12, 16, etc.
             */
            private int numConfsPerBatch = 8;

            public Builder setAllowedSeqs(EwakstarLimitedSequenceTrie seqs) {
                allowedSeqs = seqs;
                return this;
            }

            public Builder setNumConfsPerBatch(int val) {
                numConfsPerBatch = val;
                return this;
            }

            public Builder setNumBestSequences(int val){
                numBestSequences = val;
                return this;
            }


            public Settings build() {
                return new Settings(allowedSeqs, numBestSequences, numConfsPerBatch);
            }
        }

        public final EwakstarLimitedSequenceTrie allowedSeqs;
        public final int numBestSequences;
        public final int numConfsPerBatch;

        public Settings(EwakstarLimitedSequenceTrie seqs, int numBestSequences, int numConfsPerBatch) {

            this.numBestSequences = numBestSequences;
            this.numConfsPerBatch = numConfsPerBatch;
            this.allowedSeqs = seqs;
        }
    }

    public class ConfSpaceInfo {

        public final SimpleConfSpace confSpace;
        public final EWAKStar.ConfSpaceType type;
        public final String id;

        /** A ConfEnergyCalculator that computes minimized energies */
        public ConfEnergyCalculator confEcalcMinimized = null;

        /** A ConfSearch that minimizes over conformation scores from minimized tuples */
        public Function<RCs,ConfAStarTree> confTreeFactoryMinimized = null;

        /** A ConfSearch that maximizes over conformation scores from rigid tuples */
        public Function<RCs,ConfAStarTree> confTreeFactoryRigid = null;

        public File confDBFile = null;

        private BigDecimal stabilityThreshold = null;

        public ConfSpaceInfo(SimpleConfSpace confSpace, EWAKStar.ConfSpaceType type) {
            this.confSpace = confSpace;
            this.type = type;
            this.id = type.name().toLowerCase();
        }

        private void check() {
            if (confEcalcMinimized == null) {
                throw new EWAKStar.InitException(type, "confEcalcMinimized");
            }
            if (confTreeFactoryMinimized == null) {
                throw new EWAKStar.InitException(type, "confTreeFactoryMinimized");
            }
            if (confTreeFactoryRigid == null) {
                throw new EWAKStar.InitException(type, "confTreeFactoryRigid");
            }
        }

        public void setConfDBFile(String path) {
            confDBFile = new File(path);
        }
    }

    private class ConfDBs {
        public ConfDB protein = null;
        public ConfDB ligand = null;
        public ConfDB complex = null;
    }

    private static interface DBsUser {
        void use(ConfDBs confdbs);
    }

    public static enum PfuncsStatus {
        Estimating,
        Estimated,
        Blocked
    }

    private abstract class Node implements Comparable<Node> {

        public final Sequence sequence;
        public final ConfDB.DBs confDBs;
        public boolean isWTSeq = false;

        /** for comparing in the tree, higher is first */
        public double score;

        /** signals whether or not partition function values are allowed the stability threshold */
        public boolean isUnboundUnstable;

        protected Node(Sequence sequence, ConfDB.DBs confDBs) {
            if (sequence.isWildType()){
                this.isWTSeq = true;
            }
            this.sequence = sequence;
            this.confDBs = confDBs;
            this.score = 0;
        }

        @Override
        public int compareTo(Node other) {
            // negate for descending sort
            return -Double.compare(this.score, other.score);
        }

        public abstract void estimateScore();
    }

    public class MultiSequenceNode extends Node {

        public MultiSequenceNode(Sequence sequence, ConfDB.DBs confdbs) {
            super(sequence, confdbs);
        }

        public List<Node> makeChildren() {

            List<Node> children = new ArrayList<>();

            // pick the next design position
            // TODO: dynamic A*?
            List<SeqSpace.Position> positions = complex.confSpace.seqSpace.positions;
            SeqSpace.Position assignPos = positions.stream()
                    .filter((pos) -> !sequence.isAssigned(pos.resNum))
                    .findFirst()
                    .orElseThrow(() -> new IllegalStateException("no design positions left to choose"));

            // get the possible assignments
            //Set<SeqSpace.ResType> resTypes;
            Set<String> resTypes;

            if(assignPos.resTypes.size() == 1)
                resTypes = new HashSet<>(getResTypeList(assignPos.resTypes));
            else {
                resTypes = filterOnPreviousSeqs();
            }

            // for each assignment...
            for (String resType : resTypes) {

                // update the sequence with this assignment
                Sequence s = sequence.copy().set(assignPos, resType);

                if (s.isFullyAssigned()) {
                    // fully assigned, make single sequence node

                    children.add(new SingleSequenceNode(s, confDBs));

                } else if (kstarSettings.useExact && s.countMutations() == kstarSettings.numMutations) {
                    // mutation limit reached, fill unassigned positions with wild-type
                    s.fillWildType();
                    children.add(new SingleSequenceNode(s, confDBs));
                } else if (s.countMutations() == kstarSettings.numMutations){
                    s.fillWildType();
                    children.add(new SingleSequenceNode(s, confDBs));
                } else {
                    // still partial sequence, make multi-sequence node
                    children.add(new MultiSequenceNode(s,confDBs));
                }
            }

            return children;
        }

        private Set<String> getResTypeList(List<SeqSpace.ResType> resTypes){

            Set<String> resTypeList = new HashSet<>();
            for(SeqSpace.ResType r: resTypes){
                resTypeList.add(r.name);
            }

            return resTypeList;
        }

        private Set<String> filterOnPreviousSeqs(){

            String subSeq = sequence.toString();
            Set<String> resTypes;
            if (subSeq.equals("")){
                 resTypes = bbkstarSettings.allowedSeqs.getFirstPos();
            }
            else {
                resTypes = bbkstarSettings.allowedSeqs.getSeq(subSeq);
            }

            return resTypes;
        }

        @Override
        public void estimateScore() {

            // TODO: expose setting?
            // NOTE: for the correctness of the bounds, the number of confs must be the same for every node
            // meaning, it might not be sound to do epsilon-based iterative approximations here

            final int numConfs = 1000;

            BigDecimal proteinLowerBound = calcLowerBound(protein, sequence, numConfs);

            // if the first few conf upper bound scores (for the pfunc lower bound) are too high,
            // then the K* upper bound is also too high
            if (MathTools.isZero(proteinLowerBound)) {
                score = Double.POSITIVE_INFINITY;
                isUnboundUnstable = false;
                return;
            }

            BigDecimal ligandLowerBound = calcLowerBound(ligand, sequence, numConfs);

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

        private BigDecimal calcLowerBound(ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute lower bounds on pfuncs, we'll do the usual lower bound calculation,
            // but use rigid energies instead of minimized energies

            RCs rcs = sequence.makeRCs(info.confSpace);
            ConfSearch astar = info.confTreeFactoryRigid.apply(rcs);
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

        private BigDecimal calcUpperBound(ConfSpaceInfo info, Sequence sequence, int numConfs) {

            // to compute upper bounds on pfuncs,
            // we'll use the upper bound calculator in the usual way

            RCs rcs = sequence.makeRCs(info.confSpace);
            UpperBoundCalculator calc = new UpperBoundCalculator(
                    info.confTreeFactoryMinimized.apply(rcs),
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

    public class SingleSequenceNode extends Node {

        public EWAKStarPartitionFunction protein;
        public EWAKStarPartitionFunction ligand;
        public EWAKStarPartitionFunction complex;

        public SingleSequenceNode(Sequence sequence, ConfDB.DBs confDBs) {
            super(sequence, confDBs);

            // make the partition functions
            this.protein = makePfunc(proteinPfuncs, EWAKStarBBKStar.this.protein, confDBs.get(EWAKStarBBKStar.this.protein.confSpace));
            this.ligand = makePfunc(ligandPfuncs, EWAKStarBBKStar.this.ligand, confDBs.get(EWAKStarBBKStar.this.ligand.confSpace));
            this.complex = makePfunc(complexPfuncs, EWAKStarBBKStar.this.complex, confDBs.get(EWAKStarBBKStar.this.complex.confSpace));
        }

        private EWAKStarPartitionFunction makePfunc(Map<Sequence,EWAKStarPartitionFunction> pfuncCache, ConfSpaceInfo info, ConfDB confdb) {

            // filter the global sequence to this conf space
            Sequence sequence = this.sequence.filter(info.confSpace.seqSpace);

            // first check the cache
            EWAKStarPartitionFunction pfunc = pfuncCache.get(sequence);
            if (pfunc != null) {
                return pfunc;
            }

            // cache miss, need to compute the partition function

            // make the partition function
            pfunc = new EWAKStarGradientDescentPfunc(info.confEcalcMinimized);
            pfunc.setReportProgress(kstarSettings.showPfuncProgress);
            if (confdb != null) {
                EWAKStarPartitionFunction.WithConfTable.setOrThrow(pfunc, confdb.getSequence(sequence));
            }
            RCs rcs = sequence.makeRCs(info.confSpace);
            if (kstarSettings.useExternalMemory) {
                EWAKStarPartitionFunction.WithExternalMemory.setOrThrow(pfunc, true, rcs);
            }
            pfunc.init(confTrees.make(() -> info.confTreeFactoryMinimized.apply(rcs)), confTrees.make(() -> info.confTreeFactoryMinimized.apply(rcs)), rcs.getNumConformations(), kstarSettings.epsilon, kstarSettings.eW, kstarSettings.maxPFConfs, kstarSettings.printPDBs);
            pfunc.setStabilityThreshold(info.stabilityThreshold);

            // update the cache
            pfuncCache.put(sequence, pfunc);
            return pfunc;
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
                protein.compute(bbkstarSettings.numConfsPerBatch);

                // tank the sequence if the unbound protein is unstable
                if (protein.getStatus() == EWAKStarPartitionFunction.Status.Unstable) {
                    score = Double.NEGATIVE_INFINITY;
                    isUnboundUnstable = true;
                    return;
                }
            }

            if (ligand.getStatus().canContinue()) {
                ligand.compute(bbkstarSettings.numConfsPerBatch);

                // tank the sequence if the unbound ligand is unstable
                if (ligand.getStatus() == EWAKStarPartitionFunction.Status.Unstable) {
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

            // tank sequences that have no useful K* bounds, and are blocked
            if (getStatus() == PfuncsStatus.Blocked && score == Double.POSITIVE_INFINITY) {
                score = Double.NEGATIVE_INFINITY;
            }
        }

        public EWAKStarScore makeKStarScore() {
            return new EWAKStarScore(protein.makeResult(), ligand.makeResult(), complex.makeResult());
        }

        public PfuncsStatus getStatus() {

            // aggregate pfunc statuses
            if ((protein.getStatus() == EWAKStarPartitionFunction.Status.ConfLimitReached
                    || protein.getStatus() == EWAKStarPartitionFunction.Status.EpsilonReached
                    || protein.getStatus() == EWAKStarPartitionFunction.Status.EnergyReached) &&
                    (ligand.getStatus() == EWAKStarPartitionFunction.Status.ConfLimitReached
                            || ligand.getStatus() == EWAKStarPartitionFunction.Status.EpsilonReached
                            || ligand.getStatus() == EWAKStarPartitionFunction.Status.EnergyReached) &&
                    (complex.getStatus() == EWAKStarPartitionFunction.Status.ConfLimitReached
                            || complex.getStatus() == EWAKStarPartitionFunction.Status.EpsilonReached
                            || complex.getStatus() == EWAKStarPartitionFunction.Status.EnergyReached)) {
                return PfuncsStatus.Estimated;
            } else if (protein.getStatus() == EWAKStarPartitionFunction.Status.Estimating
                    || ligand.getStatus() == EWAKStarPartitionFunction.Status.Estimating
                    || complex.getStatus() == EWAKStarPartitionFunction.Status.Estimating) {
                return PfuncsStatus.Estimating;
            } else {
                return PfuncsStatus.Blocked;
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
    public final ConfSpaceInfo protein;

    /** A configuration space containing just the ligand strand */
    public final ConfSpaceInfo ligand;

    /** A configuration space containing both the protein and ligand strands */
    public final ConfSpaceInfo complex;

    /** Optional and overridable settings for BBK*, shared with K* */
    public final EWAKStar.Settings kstarSettings;
    public final Settings bbkstarSettings;

    // TODO: caching these will keep lots of A* trees in memory. is that a problem?
    private final Map<Sequence,EWAKStarPartitionFunction> proteinPfuncs;
    private final Map<Sequence,EWAKStarPartitionFunction> ligandPfuncs;
    private final Map<Sequence,EWAKStarPartitionFunction> complexPfuncs;

    private final ConfSearchCache confTrees;

    public EWAKStarBBKStar(EwakstarDoer.State P, EwakstarDoer.State L, EwakstarDoer.State PL, EWAKStar.Settings kstarSettings, EWAKStarBBKStar.Settings bbkstarSettings, Integer minNumConfTrees) {

        // BBK* doesn't work with external memory (never enough internal memory for all the priority queues)
        if (kstarSettings.useExternalMemory) {
            throw new IllegalArgumentException("BBK* is not compatible with external memory."
                    + " Please switch to regular K* with external memory, or keep using BBK* and disable external memory.");
        }

        this.protein = new ConfSpaceInfo(P.confSpace, EWAKStar.ConfSpaceType.Protein);
        this.ligand = new ConfSpaceInfo(L.confSpace, EWAKStar.ConfSpaceType.Ligand);
        this.complex = new ConfSpaceInfo(PL.confSpace, EWAKStar.ConfSpaceType.Complex);
        this.kstarSettings = kstarSettings;
        this.bbkstarSettings = bbkstarSettings;

        confTrees = new ConfSearchCache(minNumConfTrees);

        proteinPfuncs = new HashMap<>();
        ligandPfuncs = new HashMap<>();
        complexPfuncs = new HashMap<>();
    }

    public Iterable<ConfSpaceInfo> confSpaceInfos() {
        return Arrays.asList(protein, ligand, complex);
    }

    public ConfSpaceInfo getConfSpaceInfo(SimpleConfSpace confSpace) {
        if (confSpace == protein.confSpace) {
            return protein;
        } else if (confSpace == ligand.confSpace) {
            return ligand;
        } else if (confSpace == complex.confSpace) {
            return complex;
        } else {
            throw new IllegalArgumentException("conf space does not match any known by this K* instance");
        }
    }

    public List<Sequence> run() {

        protein.check();
        ligand.check();
        complex.check();

        // clear any previous state
        proteinPfuncs.clear();
        ligandPfuncs.clear();
        complexPfuncs.clear();

        List<Sequence> scoredSequences = new ArrayList<>();


        try (ConfDB.DBs confDBs = new ConfDB.DBs()
                .add(protein.confSpace, protein.confDBFile)
                .add(ligand.confSpace, ligand.confDBFile)
                .add(complex.confSpace, complex.confDBFile)
        ) {

            // start the BBK* tree with the root node
            PriorityQueue<Node> tree = new PriorityQueue<>();
            tree.add(new MultiSequenceNode(complex.confSpace.makeUnassignedSequence(), confDBs));

            // start searching the tree
            System.out.println("Computing K* scores for the best sequences to within an energy window of " + kstarSettings.eW + " kcal, with a max of " + kstarSettings.maxPFConfs + " conformations, and an epsilon of " + kstarSettings.epsilon + "...");
            kstarSettings.scoreWriters.writeHeader();

            Boolean wtSeqFound = false;
            if (kstarSettings.wtBenchmark) {
                while (!tree.isEmpty() && !wtSeqFound) {

                    // get the next node
                    Node node = tree.poll();

                    if (node.isWTSeq) {
                        wtSeqFound = true;
                    }

                    if (node instanceof SingleSequenceNode) {
                        SingleSequenceNode ssnode = (SingleSequenceNode) node;

                        // single-sequence node
                        switch (ssnode.getStatus()) {
                            case Estimated:

                                // sequence is finished, return it!
                                reportSequence(ssnode, scoredSequences);
                                ssnode.protein = null;
                                ssnode.complex = null;
                                ssnode.ligand = null;
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

                    } else if (node instanceof MultiSequenceNode) {
                        MultiSequenceNode msnode = (MultiSequenceNode) node;

                        // partial sequence, expand children
                        // TODO: parallelize the multi-sequence node scoring here?
                        for (Node child : msnode.makeChildren()) {
                            child.estimateScore();
                            if (!child.isUnboundUnstable) {
                                tree.add(child);
                            }
                        }
                    }
                }
            } else {
                while (!tree.isEmpty() && scoredSequences.size() < kstarSettings.numTopOverallSeqs) {

                    // get the next node
                    Node node = tree.poll();

                    if (node instanceof SingleSequenceNode) {
                        SingleSequenceNode ssnode = (SingleSequenceNode) node;

                        // single-sequence node
                        switch (ssnode.getStatus()) {
                            case Estimated:

                                // sequence is finished, return it!
                                reportSequence(ssnode, scoredSequences);
                                ssnode.protein = null;
                                ssnode.complex = null;
                                ssnode.ligand = null;
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

                    } else if (node instanceof MultiSequenceNode) {
                        MultiSequenceNode msnode = (MultiSequenceNode) node;

                        // partial sequence, expand children
                        // TODO: parallelize the multi-sequence node scoring here?
                        for (Node child : msnode.makeChildren()) {
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
            } else if (wtSeqFound && scoredSequences.size() == 1) {
                System.out.println("No K* scores found that are better than the wild-type sequence.");
            } else if (wtSeqFound) {
                System.out.println("Found K* score estimates for all " + scoredSequences.size() + " sequences with scores greater than that of the wild-type sequence.");
            } else if (scoredSequences.size() >= kstarSettings.numTopOverallSeqs) {
                System.out.println("Found K* score estimates for top " + kstarSettings.numTopOverallSeqs + " sequences.");
            } else
                throw new Error("EWAK* ended, but the tree isn't empty and we didn't return all of the sequences. This is a bug.");

        }
        return scoredSequences;
    }

    private void reportSequence(SingleSequenceNode ssnode, List<Sequence> scoredSequences) {

        EWAKStarScore kstarScore = ssnode.makeKStarScore();
        scoredSequences.add(ssnode.sequence);

        //System.out.println(ssnode.complex.getSConfs().size());

        if(kstarSettings.printPDBs) {
            Iterator<EnergyCalculator.EnergiedParametricMolecule> econfs = ssnode.complex.getEpMols().iterator();
            HashMap<Double, ConfSearch.ScoredConf> sconfs = ssnode.complex.getSConfs();

            // return the analysis
            ConfAnalyzer analyzer = new ConfAnalyzer(complex.confEcalcMinimized);
            ConfAnalyzer.EnsembleAnalysis analysis = analyzer.analyzeEnsemble(sconfs, econfs, 10);
            String pdbString = "pdbs";
            File pdbDir = new File(pdbString);
            if (!pdbDir.exists()) {
                pdbDir.mkdir();
            }
            String seqDir = ssnode.sequence.toString().replaceAll(" ", "_");
            File directory = new File(pdbString + "/" + seqDir);
            if (!directory.exists()) {
                directory.mkdir();
            }
            analysis.writePdbs(pdbString + "/" + seqDir + "/conf.*.pdb");
            sconfs = null;
            analyzer = null;
            econfs = null;
        }

        kstarSettings.scoreWriters.writeScore(new EWAKStarScoreWriter.ScoreInfo(
                scoredSequences.size() - 1,
                0,
                ssnode.sequence,
                complex.confSpace,
                kstarScore
        ));

        kstarScore = null;
    }
}
