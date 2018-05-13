/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;

import java.io.File;
import java.util.*;
import java.util.function.Function;

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
			buf.append(String.format("Residues           %s\n", confSpace.formatResidueNumbers()));
			buf.append(String.format("%-16s   %s\n", info.id + " Sequence", sequence.toString(Sequence.Renderer.ResTypeMutations, 5)));
			buf.append(String.format("Ensemble of %d conformations:\n", ensemble.analyses.size()));
			for (int i=0; i<ensemble.analyses.size(); i++) {
				ConfAnalyzer.ConfAnalysis analysis = ensemble.analyses.get(i);
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

	public static class ConfSpaceInfo {
		public SimpleConfSpace confSpace;
		public String id;
		public ConfEnergyCalculator confEcalc;
		public KStar.ConfSearchFactory confSearchFactory;
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

				if (info.confSpace != sequence.confSpace) {
					continue;
				}

				ConfSpaceInfo adapter = new ConfSpaceInfo();
				adapter.confSpace = info.confSpace;
				adapter.id = info.id;
				adapter.confEcalc = info.confEcalc;
				adapter.confSearchFactory = info.confSearchFactory;
				adapter.confDBFile = info.getConfDBFile();
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

				if (info.confSpace != sequence.confSpace) {
					continue;
				}

				ConfSpaceInfo adapter = new ConfSpaceInfo();
				adapter.confSpace = info.confSpace;
				adapter.id = info.id;
				adapter.confEcalc = info.confEcalcMinimized;
				adapter.confSearchFactory = info.confSearchFactoryMinimized;
				adapter.confDBFile = info.getConfDBFile();
				return adapter;
			}

			throw new NoSuchConfSpaceException(sequence);
		};
	}

	public Analysis analyze(Sequence sequence, double energyWindowSize) {

		ConfSpaceInfo info = finder.apply(sequence);

		// find the GMEC for this sequence
		ConfSearch astar = info.confSearchFactory.make(sequence.makeRCs());
		SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, info.confEcalc)
			.setConfDB(info.confDBFile)
			.build();
		Queue.FIFO<ConfSearch.EnergiedConf> econfs = gmecFinder.find(energyWindowSize);

		// return the analysis
		ConfAnalyzer analyzer = new ConfAnalyzer(info.confEcalc);
		ConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(econfs, Integer.MAX_VALUE);
		return new Analysis(info, sequence, ensemble);
	}
}
