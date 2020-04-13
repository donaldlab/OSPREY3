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

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;

import java.io.File;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * Shows information about a single sequence.
 */
public class SequenceAnalyzer {

	public class Analysis {

		public final ConfSpaceInfo info;
		public final Sequence sequence;
		public final ConfAnalyzer.EnsembleAnalysis ensemble;

		public Analysis(ConfSpaceInfo info, Sequence sequence, ConfAnalyzer.EnsembleAnalysis ensemble) {
			this.info = info;
			this.sequence = sequence;
			this.ensemble = ensemble;
		}

		public ConfSpaceIteration getConfSpace() {
			return info.confSpace;
		}

		public void writePdbs(String filePattern) {
			ensemble.writePdbs(filePattern);
		}

		public void writePdb(String path, String comment) {
			ensemble.writePdb(path, comment);
		}

		@Override
		public String toString() {

			ConfSpaceIteration confSpace = getConfSpace();
			int indexSize = 1 + (int)Math.log10(ensemble.analyses.size());

			StringBuilder buf = new StringBuilder();
			buf.append("Residues           "); // TODO: rename to something not residue-specific?
			buf.append(IntStream.range(0, confSpace.numPos())
				.mapToObj(posi -> String.format("%-5s", confSpace.name(posi)))
				.collect(Collectors.joining(" "))
			);
			buf.append(String.format("%-16s   %s\n", info.id + " Sequence", sequence.toString(Sequence.Renderer.ResTypeMutations, 5)));
			buf.append(String.format("Ensemble of %d conformations:\n", ensemble.analyses.size()));
			for (int i=0; i<ensemble.analyses.size(); i++) {
				ConfAnalyzer.ConfAnalysis analysis = ensemble.analyses.get(i);
				buf.append("\t");
				buf.append(String.format("%" + indexSize + "d/%" + indexSize + "d", i + 1, ensemble.analyses.size()));
				buf.append(String.format("     Energy: %-12.6f     Score: %-12.6f", analysis.epmol.energy, analysis.score));
				buf.append("     Rotamers: "); // TODO: rename to something not rotamer-specific?
				buf.append(IntStream.range(0, confSpace.numPos())
					.mapToObj(posi -> confSpace.confId(posi, analysis.assignments[posi]))
					.collect(Collectors.joining(" "))
				);
				buf.append("     Residue Conf IDs: ");
				buf.append(SimpleConfSpace.formatConfRCs(analysis.assignments));
				buf.append("\n");
			}
			return buf.toString();
		}
	}

	public static class ConfSpaceInfo {
		public ConfSpaceIteration confSpace;
		public String id;
		public ConfEnergyCalculator confEcalc;
		public File confDBFile;
	}

	public static class NoSuchConfSpaceException extends NoSuchElementException {

		public NoSuchConfSpaceException(Sequence sequence) {
			super("no conformation space matching sequence " + sequence);
		}
	}


	private final Function<Sequence,ConfSpaceInfo> finder;

	/**
	 * make a SequenceAnalyzer from a KStar instance
	 */
	public SequenceAnalyzer(KStar kstar) {

		finder = (sequence) -> {

			for (KStar.ConfSpaceInfo info : kstar.confSpaceInfos()) {

				if (info.confSpace.seqSpace() != sequence.seqSpace) {
					continue;
				}

				ConfSpaceInfo adapter = new ConfSpaceInfo();
				adapter.confSpace = info.confSpace;
				adapter.id = info.id;
				adapter.confEcalc = info.confEcalc;
				adapter.confDBFile = info.confDBFile;
				return adapter;
			}

			throw new NoSuchConfSpaceException(sequence);
		};
	}

	/**
	 * make a SequenceAnalyzer from a BBKStar instance
	 */
	public SequenceAnalyzer(BBKStar bbkstar) {

		finder = (sequence) -> {

			for (BBKStar.ConfSpaceInfo info : bbkstar.confSpaceInfos()) {

				if (info.confSpace.seqSpace() != sequence.seqSpace) {
					continue;
				}

				ConfSpaceInfo adapter = new ConfSpaceInfo();
				adapter.confSpace = info.confSpace;
				adapter.id = info.id;
				adapter.confEcalc = info.confEcalcMinimized;
				adapter.confDBFile = info.confDBFile;
				return adapter;
			}

			throw new NoSuchConfSpaceException(sequence);
		};
	}

	/**
	 * Analyzes the sequence by getting the ensemble purely from the conformation database.
	 * Will only calculate structures for the low-energy conformations in the conformation
	 * database, so guaranteed to be fast.
	 */
	public Analysis analyze(Sequence sequence, int numConfs) {
		return analyze(sequence, numConfs, null);
	}

	/**
	 * Analyzes the sequence with the given energy calculator.
	 * Useful for LUTE, whose KStar settings use the LUTE energy calculator, which can't return atomic models.
	 */
	public Analysis analyze(Sequence sequence, int numConfs, ConfEnergyCalculator confEcalc) {

		ConfSpaceInfo info = finder.apply(sequence);

		// iterate through conformations in the confdb in order of (weakly) increasing energy
		try (ConfDB confdb = new ConfDB(info.confSpace, info.confDBFile)) {

			Iterator<ConfSearch.EnergiedConf> econfs = confdb
				.getSequence(sequence)
				.energiedConfs(ConfDB.SortOrder.Energy)
				.iterator();

			if (confEcalc == null) {
				confEcalc = info.confEcalc;
			}

			// return the analysis
			ConfAnalyzer analyzer = new ConfAnalyzer(confEcalc);
			ConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(econfs, numConfs);
			return new Analysis(info, sequence, ensemble);
		}
	}
}
