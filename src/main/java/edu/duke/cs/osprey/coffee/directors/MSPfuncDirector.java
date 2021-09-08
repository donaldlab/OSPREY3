package edu.duke.cs.osprey.coffee.directors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.coffee.FreeEnergyCalculator;
import edu.duke.cs.osprey.coffee.NodeProcessor;
import edu.duke.cs.osprey.coffee.directions.Directions;
import edu.duke.cs.osprey.coffee.nodedb.NodeTree;
import edu.duke.cs.osprey.coffee.seqdb.SeqFreeEnergies;
import edu.duke.cs.osprey.coffee.seqdb.SeqInfo;
import edu.duke.cs.osprey.coffee.seqdb.StateZ;
import edu.duke.cs.osprey.confspace.MultiStateConfSpace;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.parallelism.ThreadTools;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;

/**
 * Pfunc director with support for tracking sequences separately when multiple sequences allowed
 */
public class MSPfuncDirector extends PfuncDirector{
    public static class Builder {
        public final MultiStateConfSpace confSpace;
        public final MultiStateConfSpace.State state;

        private Sequence seq;
        private double gWidthMax = 1.0;
        private Double gMax = null;
        private Timing timing = Timing.Efficient;
        private boolean reportProgress = false;
        protected int ensembleSize = 0;
        private File ensembleDir = null;
        private long ensembleUpdate = 30;
        private TimeUnit ensembleUpdateUnit = TimeUnit.SECONDS;
        private Integer K;
        private Integer maxSimultaneousMutations;

        public Builder(MultiStateConfSpace confSpace, MultiStateConfSpace.State state) {
            this.confSpace = confSpace;
            this.state = state;
        }

        public Builder(MultiStateConfSpace confSpace, String stateName) {
            this(confSpace, confSpace.getState(stateName));
        }

        public Builder(ConfSpace state) {
            this(
                    new MultiStateConfSpace
                            .Builder("state", state)
                            .build(),
                    "state"
            );
        }

        public MSPfuncDirector.Builder setSequence(Sequence val) {

            if (!state.isSequenced && seq != null) {
                log("WARNING: Ignoring sequence given for unsequenced state %s", state.name);
                seq = null;
            }

            seq = val;
            return this;
        }

        /**
         * Sets the largest precision desired for free energy calculations.
         * Resulting free energy values may be more precise than this maximum value.
         * If the computation is stopped early, the free energy values may be less precise than this value.
         */
        public MSPfuncDirector.Builder setGWidthMax(double val) {
            gWidthMax = val;
            return this;
        }

        /**
         * Sets the maximum free energy value that is useful.
         * The free energy calculation will be aborted early if the lower bound
         * can be proven to be above this maximum value.
         */
        public MSPfuncDirector.Builder setGMax(Double val) {
            gMax = val;
            return this;
        }

        public MSPfuncDirector.Builder setTiming(Timing val) {
            timing = val;
            return this;
        }

        public MSPfuncDirector.Builder setReportProgress(boolean val) {
            reportProgress = val;
            return this;
        }

        /**
         * Tracks the lowest-energy conformations for the best sequences
         * and periodically writes out ensemble PDB files to the specified directory.
         */
        public MSPfuncDirector.Builder setEnsembleTracking(int size, File dir) {
            ensembleSize = size;
            ensembleDir = dir;
            return this;
        }

        /**
         * Sets the minimum interval for writing the next lowest-energy ensemble PDB file.
         */
        public MSPfuncDirector.Builder setEnsembleMinUpdate(long update, TimeUnit updateUnit) {
            ensembleUpdate = update;
            ensembleUpdateUnit = updateUnit;
            return this;
        }


        public MSPfuncDirector.Builder setK(int val) {
            K = val;
            return this;
        }

        public MSPfuncDirector.Builder setMaxSimultaneousMutations(Integer val) {
            maxSimultaneousMutations = val;
            return this;
        }

        public MSPfuncDirector build() {
            return new MSPfuncDirector(confSpace, state, seq, K, maxSimultaneousMutations, gWidthMax, gMax, timing, reportProgress, ensembleSize, ensembleDir, ensembleUpdate, ensembleUpdateUnit);
        }

    }

    public final int K;
    public final Integer maxSimultaneousMutations;

    public final List<SeqFreeEnergies> bestSeqs;
    private final BestSet<SeqFreeEnergies> bestSeqsByLower;
    private final BestSet<SeqFreeEnergies> bestSeqsByUpper;
    private final File ensembleDir;

    protected MSPfuncDirector(MultiStateConfSpace confSpace, MultiStateConfSpace.State state, Sequence seq, Integer K, Integer maxSimultaneousMutations, double gWidthMax, Double gMax, Timing timing, boolean reportProgress, int ensembleSize, File ensembleDir, long ensembleUpdate, TimeUnit ensembleUpdateUnit) {
        super(confSpace, state, seq, gWidthMax, gMax, timing, reportProgress, ensembleSize, null, ensembleUpdate, ensembleUpdateUnit);
        this.K = K;
        this.maxSimultaneousMutations = maxSimultaneousMutations;
        this.ensembleDir = ensembleDir;

        bestSeqs = new ArrayList<>();
        bestSeqsByLower = new BestSet<>(K, s -> s.freeEnergies[state.index].lower);
        bestSeqsByUpper = new BestSet<>(K, s -> s.freeEnergies[state.index].upper);
    }


    @Override
    public void direct(Directions directions, NodeProcessor processor) {
        //TODO: Add way to restrict the sequence space based on input sequence

        directions.member.log("Searching for low free energy sequences in state %s", state.name);
        Stopwatch stopwatch = new Stopwatch().start();

        bestSeqs.clear();

        // init the cluster with the state state, all sequences
        directions.focus(state.index);
        var tree = new NodeTree(new RCs(state.confSpace), maxSimultaneousMutations);
        directions.setTree(state.index, tree);
        processor.initRootNode(state.index, tree);

        var gcalc = new FreeEnergyCalculator();
        var finishedSequences = new HashSet<Sequence>();
        var newlyFinishedSequences = new HashSet<Sequence>();

        long lastEnsembleNs = stopwatch.getTimeNs();

        // first, find the sequences with the best state free energies
        // TODO: add fail-safe to exit early if we somehow ran out of nodes to process?
        while (true) {

            // wait a bit before checking progress
            ThreadTools.sleep(timing.workMs(stopwatch), TimeUnit.MILLISECONDS);

            // get all the sequences while the seqdb is locked
            Map<Sequence, StateZ[]> seqs = new HashMap<>();
            synchronized (processor.seqdb) {
                for (var entry : processor.seqdb.getSequenced()) {
                    Sequence seq = entry.getKey();

                    SeqInfo seqInfo = entry.getValue();
                    seqs.put(seq, seqInfo.statezs);
                }
            }

            // compute all the free energies
            Map<Sequence, SeqFreeEnergies> seqgs = new HashMap<>();
            for (var entry : seqs.entrySet()) {
                Sequence seq = entry.getKey();
                StateZ[] statezs = entry.getValue();

                var seqg = new SeqFreeEnergies(seq, Arrays.stream(statezs)
                        .map(statez -> gcalc.calc(statez.zSumBounds))
                        .toArray(MathTools.DoubleBounds[]::new)
                );
                seqgs.put(seq, seqg);
            }

            // mark finished sequences
            finishedSequences.clear();
            for (var seqg : seqgs.values()) {
                var stateG = seqg.freeEnergies[state.index];
                var stateZ = seqs.get(seqg.seq)[state.index].zSumBounds;

                // finish the sequence if the Z upper bound is negative
                // but for fully assigned sequences only
                // partially-assigned sequences often have residuals that are relatively close to zero, but sometimes slightly negative
                // (indicates not enough seqdb precision to compute the bound correctly)
                if (MathTools.isNegative(stateZ.upper) && seqg.seq.isFullyAssigned()) {
                    finishedSequences.add(seqg.seq);
                }

                // finish the sequnce if the precision is on target
                if (stateG.size() <= gWidthMax) {
                    finishedSequences.add(seqg.seq);
                }
            }
            directions.finishSequences(state.sequencedIndex, finishedSequences, newlyFinishedSequences);

            // report the newly finished sequences
            for (var newlyFinishedSeq : newlyFinishedSequences) {
                var seqg = seqgs.get(newlyFinishedSeq);
                var g = seqg.freeEnergies[state.index];
                if (Double.isNaN(g.lower)) {
                    directions.member.log("WARNING: SeqDB lacks precision to compute free energy lower bound correctly for: [%s]", newlyFinishedSeq);
                }
                directions.member.log("Finished sequence: [%s]   state G %s   width %.6f   time %s",
                        newlyFinishedSeq, g, g.size(), stopwatch.getTime(2)
                );
            }

            // TODO: apply stability filters

            bestSeqs.clear();

            // find the best sequences by state free energy bounds
            bestSeqsByLower.clear();
            bestSeqsByUpper.clear();
            for (var seqg : seqgs.values()) {
                bestSeqsByLower.add(seqg);
                bestSeqsByUpper.add(seqg);
            }

            // find all the sequences that overlap the worst best upper bound
            var worstBestUpper = bestSeqsByUpper.maxScore();
            var overlapSeqgs = new ArrayList<SeqFreeEnergies>();
            for (var seqg : seqgs.values()) {
                if (seqg.freeEnergies[state.index].contains(worstBestUpper)) {
                    overlapSeqgs.add(seqg);
                }
            }

            // should we save ensembles now?
            if (ensembleUpdate > 0 && stopwatch.getTimeNs() >= lastEnsembleNs + ensembleUpdateUnit.toNanos(ensembleUpdate)) {
                saveEnsemble(directions, processor, bestSeqsByLower.toList());
                lastEnsembleNs = stopwatch.getTimeNs();
            }

            // are we done yet?

            // to be done, all the best sequences must have finite bounds
            if (!bestSeqsByUpper.toList().stream()
                    .allMatch(seqg -> seqg.freeEnergies[state.index].isFinite())
            ) {
                reportProgress(directions, processor, finishedSequences, stopwatch);
                continue;
            }

            // to be done, all the best sequences must be fully assigned
            if (!bestSeqsByUpper.toList().stream()
                    .allMatch(seqg -> seqg.seq.isFullyAssigned())
            ) {
                reportProgress(directions, processor, finishedSequences, stopwatch);
                continue;
            }

            // if all of overlaps are in the best K by upper, we're done
            if (overlapSeqgs.stream()
                    .allMatch(seqg -> bestSeqsByUpper.contains(seqg))
            ) {
                directions.member.log("Found best %d sequences!", K);

                // which means the best K by upper and the best K by lower are the same set
                var map = new HashMap<Sequence,SeqFreeEnergies>();
                for (var seqg : bestSeqsByLower.toList()) {
                    map.put(seqg.seq, seqg);
                }
                bestSeqs.addAll(map.values());

                break;
            }

            // if all the overlaps are at least the desired precision, we're done
            if (overlapSeqgs.stream()
                    .allMatch(seqg -> seqg.freeEnergies[state.index].size() <= gWidthMax)
            ) {
                directions.member.log("Can't strictly find best %d sequences, but all borderline candidates are at the desired precision, so no further progress is possible.", K);

                // in this case, the best sequences are the best K by upper, and the overlaps
                var map = new HashMap<Sequence,SeqFreeEnergies>();
                for (var seqg : bestSeqsByUpper.toList()) {
                    map.put(seqg.seq, seqg);
                }
                for (var seqg : overlapSeqgs) {
                    map.put(seqg.seq, seqg);
                }
                bestSeqs.addAll(map.values());

                break;
            }

            // nope, not done yet. write out the most promising sequences and keep going
            reportProgress(directions, processor, finishedSequences, stopwatch);
        }

        // leave a nice clean nodedb at the end
        processor.nodedb.clear(state.index);

        // all done
        directions.member.log("All states complete in %s", stopwatch.getTime(2));
        directions.member.log(report());
    }

    private void reportProgress(Directions directions, NodeProcessor processor, Set<Sequence> finishedSeqs, Stopwatch stopwatch) {

        // show best sequences that are still unfinished
        for (var score : bestSeqsByLower.scores()) {
            var seqgs = bestSeqsByLower.values(score);

            // find the unfinished sequences
            var unfinishedSeqs = seqgs.stream()
                    .filter(seqg -> !finishedSeqs.contains(seqg.seq))
                    .collect(Collectors.toList());

            if (!unfinishedSeqs.isEmpty()) {

                for (var seqg : unfinishedSeqs) {
                    var g = seqg.freeEnergies[state.index];
                    directions.member.log("Best unfinished sequence(s): [%s] complex G %s   width %.6f   finishedSeqs %d   nodedb %5.1f%%   time %s",
                            seqg.seq,
                            g.toString(3), g.size(),
                            finishedSeqs.size(),
                            processor.nodedb.usage()*100f,
                            stopwatch.getTime(2)
                    );
                }

                break;
            }
        }
    }

    public String report() {
        var buf = new StringBuilder();
        buf.append(String.format("Affinity: %d of %d sequences:\n", bestSeqsByLower.size(), K));
        bestSeqs.sort((a,b) -> a.freeEnergies[0].compareTo(b.freeEnergies[0]));
        for (var seqg : bestSeqs) {
            var stateG = seqg.freeEnergies[state.index];
            buf.append(String.format("\t[%s]   state G %s (w %.6f)\n",
                    seqg.seq,
                    stateG, stateG.size()
            ));
        }
        return buf.toString();
    }

    private void saveEnsemble(Directions directions, NodeProcessor processor, List<SeqFreeEnergies> seqgs) {

        // skip if ensemble saving is disabled
        if (ensembleDir == null) {
            return;
        }

        // keep only fully-assigned sequences
        seqgs = seqgs.stream()
                .filter(seqg -> seqg.seq.isFullyAssigned())
                .collect(Collectors.toList());

        var stopwatch = new Stopwatch().start();
        directions.member.log("Writing ensembles for %d sequences ...", seqgs.size());
        ensembleDir.mkdirs();

        for (var seqg : seqgs) {

            // pick a name for the ensemble file
            String seqstr = seqg.seq.toString(Sequence.Renderer.ResTypeMutations)
                    .replace(' ', '-');
            File ensembleFile = new File(ensembleDir, String.format("seq.%s.pdb", seqstr));

            // get the confs, if any
            List<int[]> bestConfs;
            synchronized (processor.seqdb) {
                bestConfs = processor.seqdb.getBestConfs(state, seqg.seq).stream()
                        .map(econf -> econf.getAssignments())
                        .collect(Collectors.toList());
            }
            if (bestConfs.isEmpty()) {
                // no structures, write an empty file
                PDBIO.writeFileEcoords(Collections.emptyList(), ensembleFile, "No structures for this sequence");
                continue;
            }

            // minimize them
            // TODO: optimize somehow? cache structures?
            var energiedCoords = processor.minimizeCoords(state.index, bestConfs);

            // write the PDB file
            String comment = String.format("Ensemble of %d conformations for:\n\t   State  %s\n\tSequence  [%s]",
                    energiedCoords.size(), state.name, seqg.seq.toString(Sequence.Renderer.AssignmentMutations)
            );
            PDBIO.writeFileEcoords(energiedCoords, ensembleFile, comment);
        }

        directions.member.log("Saved ensembles to %s in %s", ensembleDir.getAbsolutePath(), stopwatch.getTime(2));
    }
}
