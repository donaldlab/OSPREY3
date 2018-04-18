package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.*;
import edu.duke.cs.osprey.gmec.EWAKStarConfAnalyzer;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.SequenceAnalyzer;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Stopwatch;

import java.io.File;
import java.util.NoSuchElementException;
import java.util.stream.Stream;

/**
 *
 * @author lowegard (adapted from SequenceAnalyzer.java)
 */

public class EWAKStarSequenceAnalyzer {

    public class Analysis {

        public final EWAKStarSequenceAnalyzer.ConfSpaceInfo info;
        public final Sequence sequence;
        public final EWAKStarConfAnalyzer.EnsembleAnalysis ensemble;

        public Analysis(EWAKStarSequenceAnalyzer.ConfSpaceInfo info, Sequence sequence, EWAKStarConfAnalyzer.EnsembleAnalysis ensemble) {
            this.info = info;
            this.sequence = sequence;
            this.ensemble = ensemble;
        }

        public SimpleConfSpace getConfSpace() {
            return info.confSpace;
        }

        public void writePdbs(String filePattern) {
            ensemble.writePdbs(filePattern);
        }

        @Override
        public String toString() {

            SimpleConfSpace confSpace = getConfSpace();
            int indexSize = 1 + (int)Math.log10(ensemble.analyses.size());

            StringBuilder buf = new StringBuilder();
            buf.append(String.format("Ensemble of %d conformations:\n", ensemble.analyses.size()));
            for (int i=0; i<ensemble.analyses.size(); i++) {
                EWAKStarConfAnalyzer.ConfAnalysis analysis = ensemble.analyses.get(i);
                buf.append("\t");
                buf.append(String.format("%" + indexSize + "d/%" + indexSize + "d", i + 1, ensemble.analyses.size()));
                buf.append(String.format("     Energy: %-12.6f     Score: %-12.6f", analysis.epmol.energy, analysis.score));
                buf.append("     Rotamers: ");
                buf.append(confSpace.formatConfRotamers(analysis.assignments));
                buf.append("     Residue Conf IDs: ");
                buf.append(SimpleConfSpace.formatConfRCs(analysis.assignments));
                buf.append("\n");
            }
            return buf.toString();
        }
    }

    public class ConfSpaceInfo {
        public final SimpleConfSpace confSpace;
        public final ConfEnergyCalculator confEcalc;

        public EnergyMatrix emat = null;

        public ConfSpaceInfo(SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc, EnergyMatrix emat) {
            this.confSpace = confSpace;
            this.confEcalc = confEcalc;
            this.emat = emat;
        }
    }

    /** A configuration space containing just the protein strand */
    public final EWAKStarSequenceAnalyzer.ConfSpaceInfo protein;

    /** Calculates the energy for a conformation */
    public final ConfEnergyCalculator confECalc;

    /** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
    public final KStar.ConfSearchFactory confSearchFactory;

    /** declare the energy matrix */
    public final EnergyMatrix emat;


    public EWAKStarSequenceAnalyzer(SimpleConfSpace confSpace, ConfEnergyCalculator confECalc, PrecomputedMatrices precompMat, KStar.ConfSearchFactory confSearchFactory) {

        this.emat = precompMat.getEmat();
        this.protein = new EWAKStarSequenceAnalyzer.ConfSpaceInfo(confSpace, confECalc, emat);
        this.confECalc = confECalc;
        this.confSearchFactory = confSearchFactory;

    }

    public EWAKStarSequenceAnalyzer.Analysis analyze(Sequence sequence, double energyWindowSize) {

        EWAKStarSequenceAnalyzer.ConfSpaceInfo info = getConfSpaceFromSequence(sequence);

        // find the GMEC for this sequence
        ConfSearch astar = confSearchFactory.make(info.emat, sequence.makeRCs());
        SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, info.confEcalc)
                .build();
        Queue.FIFO<ConfSearch.EnergiedConf> econfs = gmecFinder.find(energyWindowSize);

        // return the analysis
        EWAKStarConfAnalyzer analyzer = new EWAKStarConfAnalyzer(info.confEcalc, info.emat);
        EWAKStarConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(econfs, Integer.MAX_VALUE);
        return new EWAKStarSequenceAnalyzer.Analysis(info, sequence, ensemble);
    }


    private EWAKStarSequenceAnalyzer.ConfSpaceInfo getConfSpaceFromSequence(Sequence sequence) {
        return Stream.of(protein)
                .filter((i) -> i.confSpace == sequence.confSpace)
                .findFirst()
                .orElseThrow(() -> new NoSuchElementException("no conformation space matching sequence " + sequence));
    }
}