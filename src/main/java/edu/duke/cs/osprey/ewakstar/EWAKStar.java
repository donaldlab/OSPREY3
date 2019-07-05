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

import edu.duke.cs.osprey.astar.conf.ConfSearchCache;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;

public class EWAKStar {

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

            private int numMutations = (int) Double.POSITIVE_INFINITY;

            private boolean useExact = false;

            /** For EWAKStar we want to calculate each score for a certain energy window
             * and not to an epsilon approximation
             */
            private double eW = 2.0;

            private int maxPFConfs = (int) Double.POSITIVE_INFINITY;

            private EWAKStarScoreWriter.Writers scoreWriters = new EWAKStarScoreWriter.Writers();

            /**
             * If true, prints out information to the console for each minimized conformation during
             * partition function approximation
             */
            private boolean showPfuncProgress = false;

            private boolean wtBenchmark = true;

            private boolean useExternalMemory = false;

            private int numTopOverallSeqs = (int) Double.POSITIVE_INFINITY;

            private boolean printPDBs = false;

            public Settings.Builder setPrintPDBs(boolean val){
                printPDBs = val;
                return this;
            }

            public Settings.Builder setNumTopOverallSeqs(int val){
                numTopOverallSeqs = val;
                return this;
            }

            public Settings.Builder setMaxNumConfs(int val){
                maxPFConfs = val;
                return this;
            }

            public Settings.Builder setWTBenchmark(boolean ans){
                wtBenchmark = ans;
                return this;
            }
            public Settings.Builder setEpsilon(double val){
                epsilon = val;
                return this;
            }

            public Settings.Builder setEw(double val){
                eW = val;
                return this;
            }

            public Settings.Builder setUseExact(boolean val){
                useExact = val;
                return this;
            }

            public Settings.Builder setUseExternalMemory(boolean val){
                useExternalMemory = val;
                return this;
            }

            public Settings.Builder setNumMutations(int val){
                numMutations = val;
                return this;
            }

            public Settings.Builder addScoreWriter(EWAKStarScoreWriter val) {
                scoreWriters.add(val);
                return this;
            }

            public Settings.Builder addScoreConsoleWriter(EWAKStarScoreWriter.Formatter val) {
                return addScoreWriter(new EWAKStarScoreWriter.ToConsole(val));
            }

            public Settings.Builder addScoreConsoleWriter() {
                return addScoreConsoleWriter(new EWAKStarScoreWriter.Formatter.SequenceKStarPfuncs());
            }

            public Settings.Builder addScoreFileWriter(File file, EWAKStarScoreWriter.Formatter val) {
                return addScoreWriter(new EWAKStarScoreWriter.ToFile(file, val));
            }

            public Settings.Builder addScoreFileWriter(File file) {
                return addScoreFileWriter(file, new EWAKStarScoreWriter.Formatter.Log());
            }

            public Settings.Builder setShowPfuncProgress(boolean val) {
                showPfuncProgress = val;
                return this;
            }

            public Settings build() {
                return new Settings(printPDBs, useExternalMemory, useExact, numMutations, epsilon, eW, maxPFConfs, numTopOverallSeqs, scoreWriters, showPfuncProgress, wtBenchmark);
            }
        }

        public boolean printPDBs;
        public final int numTopOverallSeqs;
        public final boolean useExternalMemory;
        public final int numMutations;
        public final boolean useExact;
        public final boolean wtBenchmark;
        public final int maxPFConfs;
        public final double epsilon;
        public final double eW;
        public final EWAKStarScoreWriter.Writers scoreWriters;
        public final boolean showPfuncProgress;

        public Settings(boolean printPDBs, boolean useExternalMemory, boolean useExact, int numMutations, double epsilon, double eW, int maxPFConfs, int numTopOverallSeqs, EWAKStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, boolean wtBenchmark) {

            this.printPDBs = printPDBs;
            this.numTopOverallSeqs = numTopOverallSeqs;
            this.useExternalMemory = useExternalMemory;
            this.useExact = useExact;
            this.numMutations = numMutations;
            this.wtBenchmark = wtBenchmark;
            this.maxPFConfs = maxPFConfs;
            this.epsilon = epsilon;
            this.eW = eW;
            this.scoreWriters = scoreWriters;
            this.showPfuncProgress = dumpPfuncConfs;
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

    public static class InitException extends RuntimeException {

        public InitException(ConfSpaceType type, String name) {
            super(String.format("set %s for the %s conf space info before running", name, type.name()));
        }
    }

    public class ConfSpaceInfo {

        public final SimpleConfSpace confSpace;
        public final ConfSpaceType type;
        public final String id;

        public final Map<Sequence,EWAKStarGradientDescentPfunc.Result> pfuncResults = new HashMap<>();

        public ConfEnergyCalculator confEcalc = null;
        public KStar.ConfSearchFactory confSearchFactory = null;
        public File confDBFile = null;

        public ConfSpaceInfo(SimpleConfSpace confSpace, ConfSpaceType type) {
            this.confSpace = confSpace;
            this.type = type;
            this.id = type.name().toLowerCase();
        }

        private void check() {
            if (confEcalc == null) {
                throw new InitException(type, "confEcalc");
            }
            if (confSearchFactory == null) {
                throw new InitException(type, "confSearchFactory");
            }
        }

        public void clear() {
            pfuncResults.clear();
        }

        public void setConfDBFile(String path) {
            confDBFile = new File(path);
        }
    }

    private static interface Scorer {
        KStarScore score(int sequenceNumber, PartitionFunction.Result proteinResult, PartitionFunction.Result ligandResult, PartitionFunction.Result complexResult);
    }

    /** A configuration space containing just the protein strand */
    public final ConfSpaceInfo protein;

    /** A configuration space containing just the ligand strand */
    public final ConfSpaceInfo ligand;

    /** A configuration space containing both the protein and ligand strands */
    public final ConfSpaceInfo complex;

    /** Optional and overridable settings for K* */
    public final Settings settings;

    public EWAKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, Settings settings) {
        this.settings = settings;
        this.protein = new ConfSpaceInfo(protein, ConfSpaceType.Protein);
        this.ligand = new ConfSpaceInfo(ligand, ConfSpaceType.Ligand);
        this.complex = new ConfSpaceInfo(complex, ConfSpaceType.Complex);

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

}
