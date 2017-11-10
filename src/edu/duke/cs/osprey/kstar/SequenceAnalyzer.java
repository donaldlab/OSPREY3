package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
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

/**
 * Shows information about a single sequence.
 */
public class SequenceAnalyzer {

	public class Analysis {

		public final KStar.ConfSpaceType type;
		public final KStar.Sequence complexSequence;
		public final KStar.Sequence filteredSequence;
		public final ConfAnalyzer.EnsembleAnalysis ensemble;

		public Analysis(KStar.ConfSpaceType type, KStar.Sequence complexSequence, KStar.Sequence filteredSequence, ConfAnalyzer.EnsembleAnalysis ensemble) {
			this.type = type;
			this.complexSequence = complexSequence;
			this.filteredSequence = filteredSequence;
			this.ensemble = ensemble;
		}

		public SimpleConfSpace getConfSpace() {
			return infosByType.get(type).confSpace;
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
			buf.append(String.format("%-16s   %s\n", type.name() + " Sequence", filteredSequence.toString(KStar.Sequence.makeWildType(confSpace), 5)));
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
		public final Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToThisMap;

		public EnergyMatrix emat = null;

		public ConfSpaceInfo(KStar.ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc, Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToThisMap) {
			this.type = type;
			this.confSpace = confSpace;
			this.confEcalc = confEcalc;
			this.complexToThisMap = complexToThisMap;
		}

		public void calcEmat() {
			SimplerEnergyMatrixCalculator.Builder builder = new SimplerEnergyMatrixCalculator.Builder(confEcalc);
			if (settings.energyMatrixCachePattern != null) {
				builder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase())));
			}
			emat = builder.build().calcEnergyMatrix();
		}

		public KStar.Sequence makeWildTypeSequence() {
			KStar.Sequence sequence = new KStar.Sequence(confSpace.positions.size());
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				sequence.set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
			}
			return sequence;
		}

		public KStar.Sequence filterSequence(KStar.Sequence complexSequence) {

			// don't need to filter complex sequences
			if (complexToThisMap == null) {
				assert (type == KStar.ConfSpaceType.Complex);
				return complexSequence;
			}

			KStar.Sequence filteredSequence = makeWildTypeSequence();
			for (SimpleConfSpace.Position pos : complex.confSpace.positions) {
				SimpleConfSpace.Position proteinPos = complexToThisMap.get(pos);
				if (proteinPos != null) {
					filteredSequence.set(proteinPos.index, complexSequence.get(pos.index));
				}
			}
			return filteredSequence;
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

	private final Map<KStar.ConfSpaceType,ConfSpaceInfo> infosByType;

	public SequenceAnalyzer(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, KStar.ConfEnergyCalculatorFactory confEcalcFactory, KStar.ConfSearchFactory confSearchFactory, KStar.Settings settings) {

		this.protein = new ConfSpaceInfo(KStar.ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc), complex.mapPositionsTo(protein));
		this.ligand = new ConfSpaceInfo(KStar.ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc), complex.mapPositionsTo(ligand));
		this.complex = new ConfSpaceInfo(KStar.ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc), null);
		this.ecalc = ecalc;
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

		// index conf space infos by type
		infosByType = new EnumMap<>(KStar.ConfSpaceType.class);
		infosByType.put(KStar.ConfSpaceType.Protein, this.protein);
		infosByType.put(KStar.ConfSpaceType.Ligand, this.ligand);
		infosByType.put(KStar.ConfSpaceType.Complex, this.complex);

		// calc emats if needed, or read them from the cache
		this.protein.calcEmat();
		this.ligand.calcEmat();
		this.complex.calcEmat();
	}

	public Analysis analyze(KStar.Sequence complexSequence, KStar.ConfSpaceType type, double energyWindowSize) {

		ConfSpaceInfo info = infosByType.get(type);

		// filter the sequence and get the RCs
		KStar.Sequence filteredSequence = info.filterSequence(complexSequence);
		RCs rcs = filteredSequence.makeRCs(info.confSpace);

		// find the GMEC for this sequence
		ConfSearch astar = confSearchFactory.make(info.emat, rcs);
		SimpleGMECFinder gmecFinder = new SimpleGMECFinder.Builder(astar, info.confEcalc).build();
		Queue.FIFO<ConfSearch.EnergiedConf> econfs = gmecFinder.find(energyWindowSize);

		// return the analysis
		ConfAnalyzer analyzer = new ConfAnalyzer(info.confEcalc, info.emat);
		ConfAnalyzer.EnsembleAnalysis ensemble = analyzer.analyzeEnsemble(econfs, Integer.MAX_VALUE);
		return new Analysis(type, complexSequence, filteredSequence, ensemble);
	}
}
