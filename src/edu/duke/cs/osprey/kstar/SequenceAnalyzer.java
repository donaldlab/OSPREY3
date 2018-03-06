package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;

import java.io.File;
import java.util.*;
import java.util.stream.Stream;

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
			buf.append(String.format("%-16s   %s\n", info.type.name() + " Sequence", sequence.toString(Sequence.Renderer.ResTypeMutations, 5)));
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

	public class ConfSpaceInfo {

		public final KStar.ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;

		public EnergyMatrix emat = null;

		public ConfSpaceInfo(KStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc) {
			this.type = type;
			this.confSpace = confSpace;
			this.confEcalc = confEcalc;
		}

		public void calcEmat() {
			SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(confEcalc);
			if (settings.energyMatrixCachePattern != null) {
				builder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase())));
			}
			emat = builder.build().calcEnergyMatrix();
		}

		public File getConfDBFile() {
			if (settings.confDBPattern == null) {
				return null;
			} else {
				return new File(settings.applyConfDBPattern(type.name().toLowerCase()));
			}
		}
	}

	/** A configuration space containing just the protein strand */
	public final ConfSpaceInfo protein;

	/** A configuration space containing just the ligand strand */
	public final ConfSpaceInfo ligand;

	/** A configuration space containing both the protein and ligand strands */
	public final ConfSpaceInfo complex;

	/** Calculates the energy for a molecule */
	public final EnergyCalculator ecalc;

	/** A function that makes a ConfEnergyCalculator with the desired options */
	public final KStar.ConfEnergyCalculatorFactory confEcalcFactory;

	/** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
	public final KStar.ConfSearchFactory confSearchFactory;

	/** K* settings */
	public final KStar.Settings settings;

	public SequenceAnalyzer(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, KStar.ConfEnergyCalculatorFactory confEcalcFactory, KStar.ConfSearchFactory confSearchFactory, KStar.Settings settings) {

		this.protein = new ConfSpaceInfo(KStar.ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc));
		this.ligand = new ConfSpaceInfo(KStar.ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc));
		this.complex = new ConfSpaceInfo(KStar.ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc));
		this.ecalc = ecalc;
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

		// calc emats if needed, or read them from the cache
		this.protein.calcEmat();
		this.ligand.calcEmat();
		this.complex.calcEmat();
	}

	public Analysis analyze(Sequence sequence, double energyWindowSize) {

		ConfSpaceInfo info = getConfSpaceFromSequence(sequence);

		// find the GMEC for this sequence
		ConfSearch astar = confSearchFactory.make(info.emat, sequence.makeRCs());
		SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, info.confEcalc)
			.setConfDB(info.getConfDBFile())
			.build();
		Queue.FIFO<ConfSearch.EnergiedConf> econfs = gmecFinder.find(energyWindowSize);

		// return the analysis
		ConfAnalyzer analyzer = new ConfAnalyzer(info.confEcalc, info.emat);
		ConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(econfs, Integer.MAX_VALUE);
		return new Analysis(info, sequence, ensemble);
	}

	public Analysis analyzeFromConfDB(Sequence sequence, double energyWindowSize) {

		ConfSpaceInfo info = getConfSpaceFromSequence(sequence);

		return new ConfDB(info.confSpace, info.getConfDBFile()).use((confdb) -> {

			// get the econfs from the table
			ConfDB.SequenceDB sdb = confdb.getSequence(sequence);
			if (sdb == null) {
				return null;
			}

			// collect the confs in the energy window
			Double minEnergy = null;
			Queue.FIFO<ConfSearch.EnergiedConf> queue = Queue.FIFOFactory.of();
			for (ConfSearch.EnergiedConf econf : sdb.energiedConfs(ConfDB.SortOrder.Energy)) {
				if (minEnergy == null) {
					minEnergy = econf.getEnergy();
				}
				if (econf.getEnergy() <= minEnergy + energyWindowSize) {
					queue.push(econf);
				}
			}

			// return the analysis
			ConfAnalyzer analyzer = new ConfAnalyzer(info.confEcalc, info.emat);
			ConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(queue, Integer.MAX_VALUE);
			return new Analysis(info, sequence, ensemble);
		});
	}

	private ConfSpaceInfo getConfSpaceFromSequence(Sequence sequence) {
		return Stream.of(complex, protein, ligand)
			.filter((i) -> i.confSpace == sequence.confSpace)
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("no conformation space matching sequence " + sequence));
	}
}
