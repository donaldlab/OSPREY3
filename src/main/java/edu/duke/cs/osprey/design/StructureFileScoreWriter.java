package edu.duke.cs.osprey.design;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.SequenceAnalyzer;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class StructureFileScoreWriter implements KStarScoreWriter {

    private final String saveDir;
    private final int numConfs;

    public StructureFileScoreWriter(String saveDir, int numConfs) {

        this.saveDir = saveDir;
        this.numConfs = numConfs;
    }

    @Override
    public void writeHeader() {

    }

    @Override
    public void writeScore(ScoreInfo info) {
        var variances = info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1).trim();
        var analysis = new SequenceAnalyzer(info.kstar).analyze(info.sequence, numConfs);
        var pdb = analysis.ensemble.writePdbString(String.format("Top %d confs for sequence", numConfs));
        var target = Paths.get(saveDir, String.format("%s.pdb", variances));

        try {
            if (!Files.exists(target.getParent())) {
                Files.createDirectory(target.getParent());
            }

            Files.writeString(target, pdb);
        } catch (IOException e) {
            System.err.printf("Could not write PDB to file %s%n", target);
            e.printStackTrace();
        }
    }
}
