package edu.duke.cs.osprey.kstar;

import java.io.File;
import java.time.Duration;

// *sigh* Java makes this stuff so verbose to do...
// Kotlin would make this so much easier
public class KStarSettings {

    public final int maxNumConfs;

    public static class Builder {

        /**
         * Value of epsilon in (0,1] for the epsilon-approximation to a partition function.
         * <p>
         * Smaller values for epsilon yield more accurate predictions, but can take
         * longer to run.
         */
        private double epsilon = 0.683;

        /**
         * Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence.
         * Defined in units of kcal/mol.
         * <p>
         * More precisely, a sequence is pruned when the following expression is true:
         * <p>
         * U(Z_s) < L(W_s) * B(t)
         * <p>
         * where:
         * - s represents the unbound protein strand, or unbound ligand strand
         * - U(Z_s) is the upper bound on the partition function for strand s
         * - L(W_s) is the lower bound on the partition function for strand s in the wild type
         * - t is the stability threshold
         * - B() is the Boltzmann weighting function
         * <p>
         * Set to null to disable the filter entirely.
         */
        private Double stabilityThreshold = 5.0;

        /**
         * The maximum number of simultaneous residue mutations to consider for each sequence mutant
         */
        private int maxSimultaneousMutations = 1;

        private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();

        /**
         * If true, prints out information to the console for each minimized conformation during
         * partition function approximation
         */
        private boolean showPfuncProgress = false;

        /**
         * True to use external memory when buffering conformations between the
         * partition function lower and upper bound calculators.
         */
        private boolean useExternalMemory = false;

        /**
         * Pattern for ConfDB filename, where the first %s is replaced with the state name.
         */
        private String confDBPattern = "%s.confdb";

        /**
         * True to attempt to resume a previous design using the conformation databases.
         * False to delete any existing conformation databases and start the design from scratch.
         */
        private boolean resume = false;

        /**
         * The maximum number of conformations for each pfunc to explore. Breaks provability in such that
         * you may not compute to a full epsilon, but can be used as a faster heuristic.
         */
        private int maxNumberConfs = Integer.MAX_VALUE;

        /**
         * The maximum amount of time to spend on any one partition function calculation.
         */
        private Duration pfuncTimeout = null;

        public Builder setEpsilon(double val) {
            epsilon = val;
            return this;
        }

        public Builder setStabilityThreshold(Double val) {
            if (val != null && val.isInfinite()) {
                throw new IllegalArgumentException("only finite values allowed. To turn off the filter, pass null");
            }
            stabilityThreshold = val;
            return this;
        }

        public Builder setMaxSimultaneousMutations(int val) {
            maxSimultaneousMutations = val;
            return this;
        }

        public Builder addScoreWriter(KStarScoreWriter val) {
            scoreWriters.add(val);
            return this;
        }

        public Builder addScoreConsoleWriter(KStarScoreWriter.Formatter val) {
            return addScoreWriter(new KStarScoreWriter.ToConsole(val));
        }

        public Builder addScoreConsoleWriter() {
            return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequenceKStarPfuncs());
        }

        public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
            return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
        }

        public Builder addScoreFileWriter(File file) {
            return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
        }

        public Builder setShowPfuncProgress(boolean val) {
            showPfuncProgress = val;
            return this;
        }

        public Builder setExternalMemory(boolean val) {
            useExternalMemory = val;
            return this;
        }

        public Builder setConfDBPattern(String val) {
            confDBPattern = val;
            return this;
        }

        public Builder setMaxNumConf(int val) {
            this.maxNumberConfs = val;
            return this;
        }

        public Builder resume(boolean val) {
            resume = val;
            return this;
        }

        public Builder setPfuncTimeout(Duration val) {
            pfuncTimeout = val;
            return this;
        }

        public KStarSettings build() {
            return new KStarSettings(epsilon, stabilityThreshold, maxSimultaneousMutations, scoreWriters, showPfuncProgress, useExternalMemory, confDBPattern, resume, maxNumberConfs, pfuncTimeout);
        }
    }

    public final double epsilon;
    public final Double stabilityThreshold;
    public final int maxSimultaneousMutations;
    public final KStarScoreWriter.Writers scoreWriters;
    public final boolean showPfuncProgress;
    public final boolean useExternalMemory;
    public final String confDBPattern;
    public final boolean resume;
    public final Duration pfuncTimeout;

    public KStarSettings(double epsilon, Double stabilityThreshold, int maxSimultaneousMutations, KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, boolean useExternalMemory, String confDBPattern, boolean resume, int maxNumberConfs, Duration pfuncTimeout) {
        this.epsilon = epsilon;
        this.stabilityThreshold = stabilityThreshold;
        this.maxSimultaneousMutations = maxSimultaneousMutations;
        this.scoreWriters = scoreWriters;
        this.showPfuncProgress = dumpPfuncConfs;
        this.useExternalMemory = useExternalMemory;
        this.confDBPattern = confDBPattern;
        this.resume = resume;
        this.maxNumConfs = maxNumberConfs;
        this.pfuncTimeout = pfuncTimeout;
    }
}
