package edu.duke.cs.osprey.design.io;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.kstar.ConfSpaceInfo;
import edu.duke.cs.osprey.kstar.NewKStar;
import edu.duke.cs.osprey.kstar.ScoredSequence;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Write the low-energy structural ensemble of a sequence
 */
public class ScoredSequencePFuncStatsWriter implements NewKStar.SequenceComputedListener {

    private final String outputDir;

    public ScoredSequencePFuncStatsWriter(String outputDir) {
        this.outputDir = outputDir;
    }

    @Override
    public void onSequence(NewKStar kstar, ScoredSequence seq) {

        // write the low-energy ensemble stats for complex, target, design
        for (Map.Entry<String, ConfSpaceInfo> entry :
                Map.of("complex.csv", kstar.complex, "design.csv", kstar.ligand, "protein.csv", kstar.protein).entrySet()) {

            var csvFile = entry.getKey();
            var info = entry.getValue();

            var seqStr = seq.sequence().toString(Sequence.Renderer.AssignmentMutations);
            var safeSeqStr = seqStr.replace(' ', '-');
            var outputFile = new File(outputDir, String.format("seq.%s.%s", safeSeqStr, csvFile));

            try (var pw = new PrintWriter(outputFile); var confDb = new ConfDB(info.confSpace, info.confDBFile)) {

                var headings = String.join(",", "score", "energy", "assignments");
                pw.println(headings);

                var localSeq = seq.sequence().filter(info.confSpace.seqSpace());
                for (ConfSearch.EnergiedConf conf : confDb.getSequence(localSeq)
                        .energiedConfs(ConfDB.SortOrder.Energy)) {
                    var energy = conf.getEnergy();
                    var score = conf.getScore();
                    var asgStr = String.join(":", Arrays.stream(conf.getAssignments()).mapToObj(Integer::toString).toList());
                    pw.println(String.format("%.6f,%.6f,%s", score, energy, asgStr));
                }
            } catch (FileNotFoundException e) {
                throw new RuntimeException(e);
            }
        }
    }
}
