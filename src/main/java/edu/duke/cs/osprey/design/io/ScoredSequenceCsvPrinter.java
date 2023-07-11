package edu.duke.cs.osprey.design.io;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.NewKStar;
import edu.duke.cs.osprey.kstar.ScoredSequence;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Writes the results of the scored sequences in a CSV format to a filestream.
 */
public class ScoredSequenceCsvPrinter implements NewKStar.SequenceComputedListener {

    private final PrintWriter writer;

    public ScoredSequenceCsvPrinter(PrintWriter writer) {
        this.writer = writer;
    }

    private boolean headerWritten;

    record Column(String header, Function<ScoredSequence, String> extractor) {
        String extract(ScoredSequence scoredSequence) {
            return extractor.apply(scoredSequence);
        }
    }

    private int row = 1;

    private final List<Column> Columns = Arrays.asList(
            new Column("Seq #", scoredSequence -> String.format("%d", row++)),
            new Column("Assignments", scoredSequence -> scoredSequence.sequence().toString(Sequence.Renderer.AssignmentMutations)),
            new Column("K* score (log10)", scoredSequence -> scoredSequence.score().scoreLog10String()),
            new Column("K* score lower (log10)", scoredSequence -> scoredSequence.score().lowerBoundLog10String()),
            new Column("K* score upper (log10)", scoredSequence -> scoredSequence.score().upperBoundLog10String()),
            new Column("Complex pfunc lower (log10)", scoredSequence -> scoredSequence.score().complexLowerBoundLog10String()),
            new Column("Complex pfunc upper (log10)", scoredSequence -> scoredSequence.score().complexUpperBoundLog10String()),
            new Column("Complex # Confs", scoredSequence -> String.format("%d", scoredSequence.score().complex.numConfs)),
            new Column("Complex Delta", scoredSequence -> String.format("%f", scoredSequence.score().complex.values.getEffectiveEpsilon())),
            new Column("Ligand pfunc lower (log10)", scoredSequence -> scoredSequence.score().ligandLowerBoundLog10String()),
            new Column("Ligand pfunc upper (log10)", scoredSequence -> scoredSequence.score().ligandUpperBoundLog10String()),
            new Column("Ligand # Confs", scoredSequence -> String.format("%d", scoredSequence.score().ligand.numConfs)),
            new Column("Ligand Delta", scoredSequence -> String.format("%f", scoredSequence.score().ligand.values.getEffectiveEpsilon())),
            new Column("Protein pfunc lower (log10)", scoredSequence -> scoredSequence.score().proteinLowerBoundLog10String()),
            new Column("Protein pfunc upper (log10)", scoredSequence -> scoredSequence.score().proteinUpperBoundLog10String()),
            new Column("Protein # Confs", scoredSequence -> String.format("%d", scoredSequence.score().protein.numConfs)),
            new Column("Protein Delta", scoredSequence -> String.format("%f", scoredSequence.score().protein.values.getEffectiveEpsilon()))
    );

    /**
     * Print a sequence in CSV format. The caller is responsible for opening and closing the stream.
     * @param sequence the sequence to write. If it's the first sequence, then the CSV column headers are written, too.
     */
    public void writeSequence(ScoredSequence sequence) {
        if (!headerWritten) {
            headerWritten = true;
            var header = Columns.stream().map(col -> col.header).collect(Collectors.joining(","));
            writer.write(header + '\n');
        }

        var rowStr = Columns.stream().map(col -> col.extract(sequence)).collect(Collectors.joining(","));
        writer.write(rowStr + '\n');
        writer.flush();
    }

    @Override
    public void onSequence(NewKStar kstar, ScoredSequence seq) {
        writeSequence(seq);
    }
}
