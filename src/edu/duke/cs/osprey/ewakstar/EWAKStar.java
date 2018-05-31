package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.io.File;
import java.util.*;

/**
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 */
public class EWAKStar {

    public interface ConfEnergyCalculatorFactory {
        ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc, String matrix);
    }

    public interface ConfSearchFactory {
        ConfSearch make(EnergyMatrix emat, RCs rcs);
    }

    // *sigh* Java makes this stuff so verbose to do...
    // Kotlin would make this so much easier
    public static class Settings {

        public static class Builder {

            /**
             * Value of epsilon in (0,1] for the epsilon-approximation to a partition function.
             *
             * Smaller values for epsilon yield more accurate predictions, but can take
             * longer to run.
             */
            private double epsilon = 0.683;

            private int maxSimultaneousMutations;

            private boolean useExact;

            /** For EWAKStar we want to calculate each score for a certain energy window
             * and not to an epsilon approximation
             */
            private double eW = 2.0;

            private int maxPFConfs = 1000000;

            private EWAKStarScoreWriter.Writers scoreWriters = new EWAKStarScoreWriter.Writers();

            /**
             * If true, prints out information to the console for each minimized conformation during
             * partition function approximation
             */
            private boolean showPfuncProgress = false;

            private boolean wtBenchmark = true;

            /**
             * Pattern of the filename to cache energy matrices.
             *
             * K*-type algorithms must calculate multiple energy matrices.
             * By default, these energy matrices are not cached between runs.
             * To cache energy matrices between runs, supply a pattern such as:
             *
             * "theFolder/emat.*.dat"
             *
             * The * in the pattern is a wildcard character that will be replaced with
             * each type of energy matrix used by the K*-type algorithm.
             */
            private String energyMatrixCachePattern = null;

            public EWAKStar.Settings.Builder setMaxNumConfs(int val){
                maxPFConfs = val;
                return this;
            }

            public EWAKStar.Settings.Builder setWTBenchmark(boolean ans){
                wtBenchmark = ans;
                return this;
            }
            public EWAKStar.Settings.Builder setEpsilon(double val){
                epsilon = val;
                return this;
            }

            public EWAKStar.Settings.Builder setEw(double val){
                eW = val;
                return this;
            }

            public EWAKStar.Settings.Builder setUseExact(boolean val){
                useExact = val;
                return this;
            }

            public EWAKStar.Settings.Builder setMaxMutations(int val){
                maxSimultaneousMutations = val;
                return this;
            }

            public EWAKStar.Settings.Builder addScoreWriter(EWAKStarScoreWriter val) {
                scoreWriters.add(val);
                return this;
            }

            public EWAKStar.Settings.Builder addScoreConsoleWriter(EWAKStarScoreWriter.Formatter val) {
                return addScoreWriter(new EWAKStarScoreWriter.ToConsole(val));
            }

            public EWAKStar.Settings.Builder addScoreConsoleWriter() {
                return addScoreConsoleWriter(new EWAKStarScoreWriter.Formatter.SequenceKStarPfuncs());
            }

            public EWAKStar.Settings.Builder addScoreFileWriter(File file, EWAKStarScoreWriter.Formatter val) {
                return addScoreWriter(new EWAKStarScoreWriter.ToFile(file, val));
            }

            public EWAKStar.Settings.Builder addScoreFileWriter(File file) {
                return addScoreFileWriter(file, new EWAKStarScoreWriter.Formatter.Log());
            }

            public EWAKStar.Settings.Builder setShowPfuncProgress(boolean val) {
                showPfuncProgress = val;
                return this;
            }

            public EWAKStar.Settings.Builder setEnergyMatrixCachePattern(String val) {
                energyMatrixCachePattern = val;
                return this;
            }

            public EWAKStar.Settings build() {
                return new EWAKStar.Settings(useExact, maxSimultaneousMutations, epsilon, eW, maxPFConfs, scoreWriters, showPfuncProgress, energyMatrixCachePattern, wtBenchmark);
            }
        }

        public final int maxSimultaneousMutations;
        public final boolean useExact;
        public final boolean wtBenchmark;
        public final int maxPFConfs;
        public final double epsilon;
        public final double eW;
        public final EWAKStarScoreWriter.Writers scoreWriters;
        public final boolean showPfuncProgress;
        public final String energyMatrixCachePattern;

        public Settings(boolean useExact, int maxSimultaneousMutations, double epsilon, double eW, int maxPFConfs, EWAKStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, String energyMatrixCachePattern, boolean wtBenchmark) {

            this.useExact = useExact;
            this.maxSimultaneousMutations = maxSimultaneousMutations;
            this.wtBenchmark = wtBenchmark;
            this.maxPFConfs = maxPFConfs;
            this.epsilon = epsilon;
            this.eW = eW;
            this.scoreWriters = scoreWriters;
            this.showPfuncProgress = dumpPfuncConfs;
            this.energyMatrixCachePattern = energyMatrixCachePattern;
        }

        public String applyEnergyMatrixCachePattern(String type) {

            // the pattern has a * right?
            if (energyMatrixCachePattern.indexOf('*') < 0) {
                throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + energyMatrixCachePattern + "') has no wildcard character (which is *)");
            }

            return energyMatrixCachePattern.replace("*", type);
        }
    }

    public static class ScoredSequence {

        public final Sequence sequence;
        public final EWAKStarScore score;

        public ScoredSequence(Sequence sequence, EWAKStarScore score) {
            this.sequence = sequence;
            this.score = score;
        }

        @Override
        public String toString() {
            return "sequence: " + sequence + "   K*(log10): " + score;
        }

        public String toString(Sequence wildtype) {
            return "sequence: " + sequence.toString(Sequence.Renderer.AssignmentMutations) + "   K*(log10): " + score;
        }
    }

    public static enum ConfSpaceType {
        Protein,
        Ligand,
        Complex
    }

    public class ConfSpaceInfo {

        public final EWAKStar.ConfSpaceType type;
        public final SimpleConfSpace confSpace;
        public final ConfEnergyCalculator confEcalc;

        public final List<Sequence> sequences = new ArrayList<>();
        public EnergyMatrix emat = null;
        public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

        public ConfSpaceInfo(EWAKStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc) {
            this.type = type;
            this.confSpace = confSpace;
            this.confEcalc = confEcalc;
        }

    }

    /** A configuration space containing just the protein strand */
    public final EWAKStar.ConfSpaceInfo protein;

    /** A configuration space containing just the ligand strand */
    public final EWAKStar.ConfSpaceInfo ligand;

    /** A configuration space containing both the protein and ligand strands */
    public final EWAKStar.ConfSpaceInfo complex;

    /** Calculates the energy for a molecule */
    public final EnergyCalculator ecalc;

    /** A function that makes a ConfEnergyCalculator with the desired options */
    public final EWAKStar.ConfEnergyCalculatorFactory confEcalcFactory;

    /** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
    public final EWAKStar.ConfSearchFactory confSearchFactory;

    /** Optional and overridable settings for K* */
    public final EWAKStar.Settings settings;

	public EWAKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, EWAKStar.ConfEnergyCalculatorFactory confEcalcFactory, EWAKStar.ConfSearchFactory confSearchFactory, EWAKStar.Settings settings, String matrixP, String matrixL, String matrixPL) {
        this.protein = new EWAKStar.ConfSpaceInfo(EWAKStar.ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc, matrixP));
        this.ligand = new EWAKStar.ConfSpaceInfo(EWAKStar.ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc, matrixL));
        this.complex = new EWAKStar.ConfSpaceInfo(EWAKStar.ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc, matrixPL));
        this.ecalc = ecalc;
        this.confEcalcFactory = confEcalcFactory;
        this.confSearchFactory = confSearchFactory;
        this.settings = settings;
    }
}
