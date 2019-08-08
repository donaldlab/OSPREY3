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

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.tools.TimeFormatter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public interface PasteScoreWriter {

	public static class ScoreInfo {

		public final int sequenceNumber;
		public final int numSequences;
		public final Sequence sequence;
		public final PasteScore pasteScore;
		public final long timeNs;

		public ScoreInfo(int sequenceNumber, int numSequences, Sequence sequence, PasteScore pasteScore) {
			this.sequenceNumber = sequenceNumber;
			this.numSequences = numSequences;
			this.sequence = sequence;
			this.pasteScore = pasteScore;
			this.timeNs = System.nanoTime();
		}
	}

	public void writeHeader();
	public void writeScore(ScoreInfo info);

	public static class Writers extends ArrayList<PasteScoreWriter> {

		private static final long serialVersionUID = 1239885431627352405L;

		public void writeHeader() {
			for (PasteScoreWriter writer : this) {
				writer.writeHeader();
			}
		}

		public void writeScore(ScoreInfo info) {
			for (PasteScoreWriter writer : this) {
				writer.writeScore(info);
			}
		}
	}

	public static class Nop implements PasteScoreWriter {

		@Override
		public void writeHeader() {
			// do nothing
		}

		@Override
		public void writeScore(ScoreInfo info) {
			// do nothing
		}
	}

	public static abstract class Formatted implements PasteScoreWriter {

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
				return String.format("%s   PAStE ddG: %-34s   ",
					info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1),
					info.pasteScore.toString()
					//info.pasteScore.protein.toString(),
					//info.pasteScore.protein.numConfs,
					//info.pasteScore.protein.values.getEffectiveEpsilon()
				);
			}
		}

		public static class Log implements Formatter {

			private final long startNs = System.nanoTime();

			@Override
			public String header() {
				return String.join("\t",
					"Sequence",
					"PAStE Score (ddG)",
					"PAStE Lower Bound",
					"PAStE Upper Bound",
					"Variant Bounds",
					"WT Bounds",
					"Total # Confs.",
					"Protein Partition Function",
					"Protein Epsilon",
					"Protein # Confs.",
					"Time (sec)"
				);
			}

			@Override
			public String format(ScoreInfo info) {
				return String.join("\t",
					info.sequence.toString(Sequence.Renderer.AssignmentMutations, info.sequence.calcCellSize() + 1),
					info.pasteScore.toString(),
					info.pasteScore.lowerBoundLog10String(),
					info.pasteScore.upperBoundLog10String(),
					info.pasteScore.protein.toString(),
					info.pasteScore.wt.toString(),
					Integer.toString(info.pasteScore.protein.numConfs ),
					String.format("%e", info.pasteScore.protein.values.qstar.doubleValue()),
					Double.toString(info.pasteScore.protein.values.getEffectiveEpsilon()),
					Integer.toString(info.pasteScore.protein.numConfs)
				);
			}
		}
	}
}
