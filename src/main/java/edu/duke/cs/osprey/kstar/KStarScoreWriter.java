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

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public interface KStarScoreWriter {

	public static class ScoreInfo {

		public final int sequenceNumber;
		public final int numSequences;
		public final Sequence sequence;
		public final KStarScore kstarScore;
		public final long timeNs;

		public ScoreInfo(int sequenceNumber, int numSequences, Sequence sequence, KStarScore kstarScore) {
			this.sequenceNumber = sequenceNumber;
			this.numSequences = numSequences;
			this.sequence = sequence;
			this.kstarScore = kstarScore;
			this.timeNs = System.nanoTime();
		}
	}

	public void writeHeader();
	public void writeScore(ScoreInfo info);

	public static class Writers extends ArrayList<KStarScoreWriter> {

		private static final long serialVersionUID = 1239885431627352405L;

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

		public static class SequenceKStarPfuncs implements Formatter {

			@Override
			public String format(ScoreInfo info) {
				return String.format("sequence %4d/%4d   %s   K*(log10): %-34s   protein: %-18s, numConfs: %d, epsilon: %01.3f   ligand: %-18s, numConfs: %d, epsilon: %01.3f   complex: %-18s, numConfs: %d, epsilon: %01.3f",
					info.sequenceNumber + 1,
					info.numSequences,
					info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1),
					info.kstarScore.toString(),
					info.kstarScore.protein.toString(),
					info.kstarScore.protein.numConfs,
					info.kstarScore.protein.values.getEffectiveEpsilon(),
					info.kstarScore.ligand.toString(),
					info.kstarScore.ligand.numConfs,
					info.kstarScore.ligand.values.getEffectiveEpsilon(),
					info.kstarScore.complex.toString(),
					info.kstarScore.complex.numConfs,
					info.kstarScore.complex.values.getEffectiveEpsilon()
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
					"K* Lower Bound",
					"K* Upper Bound",
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
					info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1),
					info.kstarScore.scoreLog10String(),
					info.kstarScore.lowerBoundLog10String(),
					info.kstarScore.upperBoundLog10String(),
					Integer.toString(info.kstarScore.protein.numConfs + info.kstarScore.ligand.numConfs + info.kstarScore.complex.numConfs),
					String.format("%e", info.kstarScore.complex.values.qstar.doubleValue()),
					Double.toString(info.kstarScore.complex.values.getEffectiveEpsilon()),
					Integer.toString(info.kstarScore.complex.numConfs),
					String.format("%e", info.kstarScore.protein.values.qstar.doubleValue()),
					Double.toString(info.kstarScore.protein.values.getEffectiveEpsilon()),
					Integer.toString(info.kstarScore.protein.numConfs),
					String.format("%e", info.kstarScore.ligand.values.qstar.doubleValue()),
					Double.toString(info.kstarScore.ligand.values.getEffectiveEpsilon()),
					Integer.toString(info.kstarScore.ligand.numConfs),
					Long.toString((info.timeNs - startNs)/TimeFormatter.NSpS)
				);
			}
		}
	}
}
