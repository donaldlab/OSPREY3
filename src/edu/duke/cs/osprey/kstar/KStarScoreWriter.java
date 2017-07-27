package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.ArrayList;

public interface KStarScoreWriter {

	public static class ScoreInfo {

		public final int sequenceNumber;
		public final int numSequences;
		public final KStar.Sequence sequence;
		public final SimpleConfSpace complexConfSpace;
		public final PartitionFunction.Result proteinResult;
		public final PartitionFunction.Result ligandResult;
		public final PartitionFunction.Result complexResult;
		public final BigDecimal kstarScore;
		public final long timeNs;

		public ScoreInfo(int sequenceNumber, int numSequences, KStar.Sequence sequence, SimpleConfSpace complexConfSpace, PartitionFunction.Result proteinResult, PartitionFunction.Result ligandResult, PartitionFunction.Result complexResult, BigDecimal kstarScore) {
			this.sequenceNumber = sequenceNumber;
			this.numSequences = numSequences;
			this.sequence = sequence;
			this.complexConfSpace = complexConfSpace;
			this.proteinResult = proteinResult;
			this.ligandResult = ligandResult;
			this.complexResult = complexResult;
			this.kstarScore = kstarScore;
			this.timeNs = System.nanoTime();
		}

		public String kstarLog10() {
			if (kstarScore != null) {
				return String.format("%f", Math.log10(kstarScore.doubleValue()));
			}
			return "none";
		}
	}

	public void writeHeader();
	public void writeScore(ScoreInfo info);

	public static class Writers extends ArrayList<KStarScoreWriter> {

		public void writeHeader() {
			for (KStarScoreWriter writer : this) {
				writer.writeHeader();
			}
		}

		public void writeScore(ScoreInfo info) {
			for (KStarScoreWriter writer : this) {
				writer.writeScore(info);
			}
		}
	}

	public static class Nop implements KStarScoreWriter {

		@Override
		public void writeHeader() {
			// do nothing
		}

		@Override
		public void writeScore(ScoreInfo info) {
			// do nothing
		}
	}

	public static abstract class Formatted implements KStarScoreWriter {

		public final Formatter formatter;

		protected Formatted(Formatter formatter) {
			this.formatter = formatter;
		}

		@Override
		public void writeHeader() {
			String header = formatter.header();
			if (header != null) {
				write(header);
			}
		}

		@Override
		public void writeScore(ScoreInfo info) {
			write(formatter.format(info));
		}

		protected abstract void write(String line);
	}

	public static class ToFile extends Formatted {

		public final File file;

		private boolean started = false;

		public ToFile(File file, Formatter formatter) {
			super(formatter);
			this.file = file;
		}

		@Override
		protected void write(String line) {

			// should we start a new file or append?
			boolean append = true;
			if (!started) {
				started = true;
				append = false;
			}

			try (FileWriter out = new FileWriter(file, append)) {
				out.write(line);
				out.write("\n");
			} catch (IOException ex) {
				System.err.println("writing to file failed: " + file);
				ex.printStackTrace(System.err);
				System.err.println(line);
			}
		}
	}

	public static class ToConsole extends Formatted {

		public ToConsole(Formatter formatter) {
			super(formatter);
		}

		@Override
		protected void write(String line) {
			System.out.println(line);
		}
	}

	public static interface Formatter {

		public default String header() {
			// no header by default
			return null;
		}

		public String format(ScoreInfo info);

		public static class SequencePfuncsScore implements Formatter {

			@Override
			public String format(ScoreInfo info) {
				return String.format("sequence %4d/%4d   %s   protein: %-18s   ligand: %-18s   complex: %-18s   K*(log10): %s",
					info.sequenceNumber + 1,
					info.numSequences,
					info.sequence,
					info.proteinResult.toString(),
					info.ligandResult.toString(),
					info.complexResult.toString(),
					info.kstarLog10()
				);
			}
		}

		public static class Log implements Formatter {

			private final long startNs = System.nanoTime();

			@Override
			public String header() {
				return String.join("\t",
					"Seq ID",
					"Sequence",
					"K* Score (Log10)",
					"Total # Confs.",
					"Complex Partition Function",
					"Complex Epsilon",
					"Complex # Confs.",
					"Protein Partition Function",
					"Protein Epsilon",
					"Protein # Confs.",
					"Ligand Partition Function",
					"Ligand Epsilon",
					"Ligand # Confs.",
					"Time (sec)"
				);
			}

			@Override
			public String format(ScoreInfo info) {
				return String.join("\t",
					Integer.toString(info.sequenceNumber),
					info.sequence.toString(info.complexConfSpace.positions),
					info.kstarLog10(),
					Integer.toString(info.proteinResult.numConfs + info.ligandResult.numConfs + info.complexResult.numConfs),
					String.format("%e", info.complexResult.values.qstar.doubleValue()),
					Double.toString(info.complexResult.values.getEffectiveEpsilon()),
					Integer.toString(info.complexResult.numConfs),
					String.format("%e", info.proteinResult.values.qstar.doubleValue()),
					Double.toString(info.proteinResult.values.getEffectiveEpsilon()),
					Integer.toString(info.proteinResult.numConfs),
					String.format("%e", info.ligandResult.values.qstar.doubleValue()),
					Double.toString(info.ligandResult.values.getEffectiveEpsilon()),
					Integer.toString(info.ligandResult.numConfs),
					Long.toString((info.timeNs - startNs)/TimeFormatter.NSpS)
				);
			}
		}

		public static class Test implements Formatter {

			@Override
			public String format(ScoreInfo info) {

				// for example:
				// assertSequence(kstar,   0, "PHE ASP GLU THR PHE LYS ILE THR", 4.050547e+04, 3.766903e+30, 3.929472e+50, epsilon); // K*(log10) = 15.410836
				return String.format("assertSequence(kstar, %3d, \"%s\", %-12s, %-12s, %-12s, epsilon); // K*(log10) = %s",
					info.sequenceNumber,
					info.sequence,
					info.proteinResult.status == PartitionFunction.Status.Estimated ? String.format("%12e", info.proteinResult.values.qstar.doubleValue()) : "null",
					info.ligandResult.status == PartitionFunction.Status.Estimated ? String.format("%12e", info.ligandResult.values.qstar.doubleValue()) : "null",
					info.complexResult.status == PartitionFunction.Status.Estimated ? String.format("%12e", info.complexResult.values.qstar.doubleValue()) : "null",
					PartitionFunction.scoreToLog10String(PartitionFunction.calcKStarScore(info.proteinResult, info.ligandResult, info.complexResult))
				);
			}
		}
	}
}
