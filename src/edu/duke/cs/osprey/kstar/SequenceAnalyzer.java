package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.gmec.SimpleGMECFinder;

import java.io.File;
import java.util.*;

public class SequenceAnalyzer {

	public class Analysis {

		public final KStar.ConfSpaceType type;
		public final KStar.Sequence complexSequence;
		public final KStar.Sequence filteredSequence;
		public final List<ConfSearch.EnergiedConf> econfs;

		public Analysis(KStar.ConfSpaceType type, KStar.Sequence complexSequence, KStar.Sequence filteredSequence) {
			this.type = type;
			this.complexSequence = complexSequence;
			this.filteredSequence = filteredSequence;
			this.econfs = new ArrayList<>();
		}

		public SimpleConfSpace getConfSpace() {
			return infosByType.get(type).confSpace;
		}

		public void writePdbs(String filePattern) {

			if (econfs.isEmpty()) {
				return;
			}

			// the pattern has a * right?
			if (filePattern.indexOf('*') < 0) {
				throw new IllegalArgumentException("filePattern (which is '" + filePattern + "') has no wildcard character (which is *)");
			}

			// mkdirs
			new File(filePattern).getParentFile().mkdirs();

			ConfSpaceInfo info = infosByType.get(type);

			String indexFormat = getIndexFormat();
			for (int i=0; i<econfs.size(); i++) {

				String filename = filePattern.replace("*", String.format(indexFormat, i + 1));
				RCTuple frag = new RCTuple(econfs.get(i).getAssignments());

				info.confEcalc.tasks.submit(() -> {
					info.confEcalc.writeMinimizedStruct(frag, filename);
					return null;
				}, (ignore) -> {});
			}

			info.confEcalc.tasks.waitForFinish();
		}

		private String getIndexFormat() {
			int indexSize = 1 + (int)Math.log10(econfs.size());
			return "%0" + indexSize + "d";
		}

		@Override
		public String toString() {

			SimpleConfSpace confSpace = getConfSpace();
			String indexFormat = getIndexFormat();

			StringBuilder buf = new StringBuilder();
			buf.append(complexSequence.toString(KStar.Sequence.makeWildType(complex.confSpace)));
			buf.append("     ");
			buf.append(type.name());
			if (!filteredSequence.equals(complexSequence)) {
				buf.append("     ");
				buf.append(filteredSequence.toString(KStar.Sequence.makeWildType(confSpace)));
			}
			buf.append("\n");
			for (int i=0; i<econfs.size(); i++) {
				ConfSearch.EnergiedConf econf = econfs.get(i);
				buf.append("\t");
				buf.append(String.format(indexFormat + "/" + indexFormat, i + 1, econfs.size()));
				buf.append(String.format("     Energy: %12.6f     Score: %12.6f", econf.getEnergy(), econf.getScore()));
				buf.append("     Rotamers: ");
				for (SimpleConfSpace.Position pos : confSpace.positions) {
					SimpleConfSpace.ResidueConf resConf = pos.resConfs.get(econf.getAssignments()[pos.index]);
					buf.append(String.format(" %3s", resConf.getRotamerCode()));
				}
				buf.append("     Residue Conf IDs: ");
				for (int rc : econf.getAssignments()) {
					buf.append(String.format(" %3d", rc));
				}
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
		Analysis analysis = new Analysis(type, complexSequence, filteredSequence);
		while (!econfs.isEmpty()) {
			analysis.econfs.add(econfs.poll());
		}
		return analysis;
	}
}
