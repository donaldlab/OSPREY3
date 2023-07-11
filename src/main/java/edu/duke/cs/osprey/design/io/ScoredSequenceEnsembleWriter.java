package edu.duke.cs.osprey.design.io;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.NewKStar;
import edu.duke.cs.osprey.kstar.ScoredSequence;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.ArrayList;

/**
 * Write the low-energy structural ensemble of a sequence
 */
public class ScoredSequenceEnsembleWriter implements NewKStar.SequenceComputedListener {

    private final int ensembleSize;
    private final String outputDir;
    private final SimpleReferenceEnergies refEnergiesComplex;
    private final ConfEnergyCalculator ecalc;
    private final ConfSpace confSpaceComplex;

    public ScoredSequenceEnsembleWriter(int ensembleSize,
                                        String outputDir,
                                        SimpleReferenceEnergies refEnergiesComplex,
                                        ConfEnergyCalculator ecalc,
                                        ConfSpace confSpaceComplex
    ) {

        this.ensembleSize = ensembleSize;
        this.outputDir = outputDir;
        this.refEnergiesComplex = refEnergiesComplex;
        this.ecalc = ecalc;
        this.confSpaceComplex = confSpaceComplex;
    }

    @Override
    public void onSequence(NewKStar kstar, ScoredSequence seq) {
        try (var confDb = new ConfDB(kstar.complex.confSpace, kstar.complex.confDBFile)) {
            var iterator = confDb.getSequence(seq.sequence())
                    .energiedConfs(ConfDB.SortOrder.Energy)
                    .iterator();

            int i = 0;
            var coords = new ArrayList<ConfEnergyCalculator.EnergiedCoords>();

            while (i < ensembleSize && iterator.hasNext()) {
                var energiedConf = iterator.next();
                var assignments = energiedConf.getAssignments();
                var posInterGen = new PosInterGen(PosInterDist.DesmetEtAl1992, refEnergiesComplex);
                var energiedCoords = ecalc.minimize(assignments, posInterGen.all(confSpaceComplex, assignments));
                coords.add(energiedCoords);
                i++;
            }

            // write the PDB file
            var seqStr = seq.sequence().toString(Sequence.Renderer.AssignmentMutations);
            var comment = String.format("Ensemble of %d conformations for:\n\t   State  %s\n\tSequence  [%s]",
                    coords.size(), confSpaceComplex.name,seqStr
            );

            // pick a name for the ensemble file
            var safeSeqStr = seqStr.replace(' ', '-');
            var ensembleFile = new File(outputDir, String.format("seq.%s.pdb", safeSeqStr));
            PDBIO.writeFileEcoords(coords, ensembleFile, comment);
        }
    }
}
