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

package edu.duke.cs.osprey.paste;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.ewakstar.EWAKStarGradientDescentPfunc;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.JvmMem;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;


/**
 * Implementation of the K* algorithm to predict protein sequence mutations that improve
 * binding affinity by computing provably accurate Boltzmann-weighted ensembles
 * {@cite Lilien2008 Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson, and Bruce R. Donald, 2005.
 * A Novel Ensemble-Based Scoring and Search Algorithm for Protein Redesign and Its Application
 * to Modify the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine Adenylation Enzyme
 * In Journal of Computational Biology (vol 12. num. 6 pp. 740–761).}.
 */
public class Paste {

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

            /**
             * can also terminate when you reach a certain energy window difference
             */

            private double eW = 10.0;

            /**
             * or when you reach a certain number of conformations
             */

            private int maxNumPfConfs = 5000;

            /**
             * decide if you want to stop calculating a partition function when the upper/lower bounds no longer
             * overlap with the WT upper/lower bounds
             */

            private boolean useWindowCriterion = true;

            /**
             * read in desired mutated sequences from a file
             */

            private File mutFile = null;

            /**
             * Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence.
             * Defined in units of kcal/mol.
             *
             * More precisely, a sequence is pruned when the following expression is true:
             *
             * U(Z_s) < L(W_s) * B(t)
             *
             * where:
             *   - s represents the unbound protein strand, or unbound ligand strand
             *   - U(Z_s) is the upper bound on the partition function for strand s
             *   - L(W_s) is the lower bound on the partition function for strand s in the wild type
             *   - t is the stability threshold
             *   - B() is the Boltzmann weighting function
             *
             * Set to null to disable the filter entirely.
             */
            private Double stabilityThreshold = null;

            /** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
            private int maxSimultaneousMutations = 1;

            private PasteScoreWriter.Writers scoreWriters = new PasteScoreWriter.Writers();

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

            public Builder setPfConfs (int val){
                maxNumPfConfs = val;
                return this;
            }

            public Builder setUseWindowCriterion(boolean val){
                useWindowCriterion = val;
                return this;
            }

            public Builder setEnergy (double val){
                eW = val;
                return this;
            }

            public Builder setNumMaxPfConfs (int val){
                maxNumPfConfs = val;
                return this;
            }

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

            public Builder addScoreWriter(PasteScoreWriter val) {
                scoreWriters.add(val);
                return this;
            }

            public Builder addScoreConsoleWriter(PasteScoreWriter.Formatter val) {
                return addScoreWriter(new PasteScoreWriter.ToConsole(val));
            }

            public Builder addScoreConsoleWriter() {
                return addScoreConsoleWriter(new PasteScoreWriter.Formatter.SequenceKStarPfuncs());
            }

            public Builder addScoreFileWriter(File file, PasteScoreWriter.Formatter val) {
                return addScoreWriter(new PasteScoreWriter.ToFile(file, val));
            }

            public Builder addScoreFileWriter(File file) {
                return addScoreFileWriter(file, new PasteScoreWriter.Formatter.Log());
            }


            public Builder addMutFile(File file){
                mutFile = file;
                return this;
            }

            public Builder setShowPfuncProgress(boolean val) {
                showPfuncProgress = val;
                return this;
            }

            public Builder setExternalMemory(boolean val) {
                useExternalMemory = val;
                return this;
            }

            public Settings build() {
                return new Settings(epsilon, eW, maxNumPfConfs, stabilityThreshold, maxSimultaneousMutations, scoreWriters, showPfuncProgress, useExternalMemory, useWindowCriterion, mutFile);
            }
        }

        public final File mutFile;
        public final double epsilon;
        public final Double stabilityThreshold;
        public final int maxSimultaneousMutations;
        public final PasteScoreWriter.Writers scoreWriters;
        public final boolean showPfuncProgress;
        public final boolean useExternalMemory;
        public final double eW;
        public final int maxNumPfConfs;
        public final boolean useWindowCriterion;


        public Settings(double epsilon, double eW, int maxNumPfConfs, Double stabilityThreshold, int maxSimultaneousMutations, PasteScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, boolean useExternalMemory, boolean useWindowCriterion, File mutFile) {
            this.eW = eW;
            this.maxNumPfConfs = maxNumPfConfs;
            this.epsilon = epsilon;
            this.stabilityThreshold = stabilityThreshold;
            this.maxSimultaneousMutations = maxSimultaneousMutations;
            this.scoreWriters = scoreWriters;
            this.showPfuncProgress = dumpPfuncConfs;
            this.useExternalMemory = useExternalMemory;
            this.useWindowCriterion = useWindowCriterion;
            this.mutFile = mutFile;
        }
    }

    public static class ScoredSequence {

        public final Sequence sequence;
        public final PasteScore score;

        public ScoredSequence(Sequence sequence, PasteScore score) {
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

    public static interface ConfSearchFactory {
        public ConfSearch make(RCs rcs);
    }

    public class ConfSpaceInfo {

        public final SimpleConfSpace confSpace;
        public final ConfSpaceType type;
        public final String id;

        public ConfEnergyCalculator confEcalc = null;
        public ConfSearchFactory confSearchFactory = null;
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

        public PastePartitionFunction.Result calcPfunc(int sequenceIndex, BigDecimal stabilityThreshold, PastePartitionFunction.Result wtResult) {

            Sequence sequence = sequences.get(sequenceIndex).filter(confSpace.seqSpace);

            // cache miss, need to compute the partition function

            // make the partition function
            PastePartitionFunction pfunc = new PasteGradientDescentPfunc(confEcalc);
            pfunc.setReportProgress(settings.showPfuncProgress);
            RCs rcs = sequence.makeRCs(confSpace);
            if (settings.useExternalMemory) {
                PastePartitionFunction.WithExternalMemory.setOrThrow(pfunc, true, rcs);
            }
            ConfSearch astar = confSearchFactory.make(rcs);
            ConfSearch astar2 = confSearchFactory.make(rcs);
            pfunc.init(astar, astar2, rcs.getNumConformations(), settings.epsilon, settings.eW, wtResult, settings.useWindowCriterion);
            pfunc.setStabilityThreshold(stabilityThreshold);

            // compute it
            pfunc.compute(settings.maxNumPfConfs);

            // save the result
            PastePartitionFunction.Result result = pfunc.makeResult();

            //pfuncResults.put(sequence, result);

			/* HACKHACK: we're done using the A* tree, pfunc, etc
				and normally the garbage collector will clean them up,
				along with their off-heap resources (e.g. TPIE data structures).
				Except the garbage collector might not do it right away.
				If we try to allocate more off-heap resources before these get cleaned up,
				we might run out. So poke the garbage collector now and try to get
				it to clean up the off-heap resources right away.
			*/

            Runtime.getRuntime().gc();

            return result;
        }

        public void setConfDBFile(String path) {
            confDBFile = new File(path);
        }
    }

    private static interface Scorer {
        PasteScore score(int sequenceNumber, PastePartitionFunction.Result proteinResult, PastePartitionFunction.Result wtResult);
    }

    /** A configuration space containing just the protein strand */
    public final ConfSpaceInfo protein;

    /** Optional and overridable settings for K* */
    public final Settings settings;

    private List<Sequence> sequences;

    public Paste(SimpleConfSpace protein, Settings settings) {
        this.settings = settings;
        this.protein = new ConfSpaceInfo(protein, ConfSpaceType.Protein);
        this.sequences = new ArrayList<>();
    }

    public Iterable<ConfSpaceInfo> confSpaceInfos() {
        return Arrays.asList(protein);
    }

    public void run() {

        // check the conf space infos to make sure we have all the inputs
        protein.check();

        // reset any previous state
        sequences.clear();

        List<ScoredSequence> scores = new ArrayList<>();

        // collect all the sequences explicitly
        sequences.add(protein.confSpace.seqSpace.makeWildTypeSequence());
        if(settings.mutFile==null) {
            sequences.addAll(protein.confSpace.seqSpace.getMutants(settings.maxSimultaneousMutations, true));
        } else {
            sequences.addAll(protein.confSpace.seqSpace.getMutants(settings.mutFile));
        }

        // TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

        // now we know how many sequences there are in total
        int n = sequences.size();

        // make the sequence scorer and reporter
        Scorer scorer = (sequenceNumber, complexResult, wtResult) -> {

            PasteScore pasteScore;
            // compute the ddG PAStE score
            if (sequenceNumber == 0){
                pasteScore = null;
            } else
                pasteScore = new PasteScore(complexResult, wtResult);
            Sequence sequence = sequences.get(sequenceNumber);
            //scores.add(new ScoredSequence(sequence, pasteScore));

            // report scores
            settings.scoreWriters.writeScore(new PasteScoreWriter.ScoreInfo(
                    sequenceNumber,
                    n,
                    sequence,
                    pasteScore
            ));

            if(pasteScore.stability.equals("Mutation Increases Stability") || pasteScore.stability.equals("Affect on Stability Unclear")) {
                Iterator<EnergyCalculator.EnergiedParametricMolecule> econfs = complexResult.epMols.iterator();
                HashMap<Double, ConfSearch.ScoredConf> sconfs = complexResult.sConfs;

                // return the analysis
                ConfAnalyzer analyzer = new ConfAnalyzer(protein.confEcalc);
                ConfAnalyzer.EnsembleAnalysis analysis = analyzer.analyzeEnsemble(sconfs, econfs, 10);
                String pdbString = "pdbs";
                File pdbDir = new File(pdbString);
                if (!pdbDir.exists()) {
                    pdbDir.mkdir();
                }
                String seqDir = sequences.get(sequenceNumber).toString().replaceAll(" ", "_");
                File directory = new File(pdbString + "/" + seqDir);
                if (!directory.exists()) {
                    directory.mkdir();
                }
                analysis.writePdbs(pdbString + "/" + seqDir + "/conf.*.pdb");
            }

            return pasteScore;
        };

        System.out.println("computing PAStE scores for " + (sequences.size()-1) + " sequence(s) to epsilon = " + settings.epsilon + " ...");
        settings.scoreWriters.writeHeader();
        // TODO: progress bar?


        // compute wild type partition functions first (always at pos 0)
        PastePartitionFunction.Result wtResult = protein.calcPfunc(0, BigDecimal.ZERO, null);
        wtResult.clearSomeResults();

        // compute all the partition functions and K* scores for the rest of the sequences
        for (int i=1; i<n; i++) {

            // get the pfuncs, with short circuits as needed
            final PastePartitionFunction.Result proteinResult = protein.calcPfunc(i, null, wtResult);

            scorer.score(i, proteinResult, wtResult);
            proteinResult.clearSomeResults();

        }
    }
}
