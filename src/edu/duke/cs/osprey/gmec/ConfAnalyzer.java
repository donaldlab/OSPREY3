package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueForcefieldBreakdown;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.structure.PDBIO;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Shows information about a single conformation.
 */
public class ConfAnalyzer {

	public class ConfAnalysis {

		public final int[] assignments;
		public final double score;
		public final EnergyCalculator.EnergiedParametricMolecule epmol;

		public ConfAnalysis(int[] assignments, double score, EnergyCalculator.EnergiedParametricMolecule epmol) {
			this.assignments = assignments;
			this.score = score;
			this.epmol = epmol;
		}

		public ResidueForcefieldBreakdown.ByResidue breakdownEnergyByResidue() {
			return new ResidueForcefieldBreakdown.ByResidue(confEcalc, epmol);
		}

		public EnergyMatrix breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type type) {
			return new ResidueForcefieldBreakdown.ByPosition(confEcalc, assignments, epmol).breakdownForcefield(type);
		}

		public EnergyMatrix breakdownScoreByPosition() {
			return new ResidueForcefieldBreakdown.ByPosition(confEcalc, assignments, epmol).breakdownScore(emat);
		}

		@Override
		public String toString() {
			return String.format(
				  "Residues           %s\n"
				+ "Sequence           %s\n"
				+ "Rotamers           %s\n"
				+ "Residue Conf IDs   %s\n"
				+ "Energy             %-12.6f\n"
				+ "Score              %-12.6f (gap: %.6f)\n",
				confEcalc.confSpace.formatResidueNumbers(),
				confEcalc.confSpace.formatConfSequence(assignments),
				confEcalc.confSpace.formatConfRotamers(assignments),
				SimpleConfSpace.formatConfRCs(assignments),
				epmol.energy,
				score,
				epmol.energy - score
			);
		}
	}

	public class EnsembleAnalysis {

		public final List<ConfAnalysis> analyses = new ArrayList<>();

		public void writePdbs(String filePattern) {
			List<EnergyCalculator.EnergiedParametricMolecule> epmols = analyses.stream()
				.map((analysis) -> analysis.epmol)
				.collect(Collectors.toList());
			PDBIO.writeEnsemble(epmols, filePattern);
		}

		@Override
		public String toString() {

			int indexSize = 1 + (int)Math.log10(analyses.size());

			StringBuilder buf = new StringBuilder();
			buf.append(String.format("Ensemble of %d conformations:\n", analyses.size()));
			buf.append("\tResidues: ");
			buf.append(confEcalc.confSpace.formatResidueNumbers());
			buf.append("\n");
			for (int i=0; i<analyses.size(); i++) {
				ConfAnalyzer.ConfAnalysis analysis = analyses.get(i);
				buf.append("\t");
				buf.append(String.format("%" + indexSize + "d/%" + indexSize + "d", i + 1, analyses.size()));
				buf.append(String.format("     Energy: %-12.6f     Score: %-12.6f", analysis.epmol.energy, analysis.score));
				buf.append("     Sequence: ");
				buf.append(confEcalc.confSpace.formatConfSequence(analysis.assignments));
				buf.append("     Rotamers: ");
				buf.append(confEcalc.confSpace.formatConfRotamers(analysis.assignments));
				buf.append("     Residue Conf IDs: ");
				buf.append(SimpleConfSpace.formatConfRCs(analysis.assignments));
				buf.append("\n");
			}
			return buf.toString();
		}
	}


	public final ConfEnergyCalculator confEcalc;
	public final EnergyMatrix emat;

	private final PairwiseGScorer gscorer;

	public ConfAnalyzer(ConfEnergyCalculator confEcalc, EnergyMatrix emat) {

		this.confEcalc = confEcalc;
		this.emat = emat;

		this.gscorer = new PairwiseGScorer(emat);
	}

	public ConfAnalysis analyze(ConfSearch.ScoredConf conf) {
		return new ConfAnalysis(
			conf.getAssignments(),
			conf.getScore(),
			confEcalc.calcEnergy(new RCTuple(conf.getAssignments()))
		);
	}

	public ConfAnalysis analyze(int[] assignments) {
		return new ConfAnalysis(
			assignments,
			gscorer.calc(assignments),
			confEcalc.calcEnergy(new RCTuple(assignments))
		);
	}

	public EnsembleAnalysis analyzeGMECEnsembleFromConfDB(String confDBPath, int maxNumConfs) {
		return analyzeGMECEnsembleFromConfDB(new File(confDBPath), maxNumConfs);
	}

	public EnsembleAnalysis analyzeGMECEnsembleFromConfDB(File confDBFile, int maxNumConfs) {
		return analyzeEnsembleFromConfDB(confDBFile, SimpleGMECFinder.ConfDBTableName, maxNumConfs);
	}

	public EnsembleAnalysis analyzeEnsembleFromConfDB(File confDBFile, String tableName, int maxNumConfs) {

		return new ConfDB(confEcalc.confSpace, confDBFile).use((confdb) -> {
			ConfDB.ConfTable table = confdb.new ConfTable(tableName);

			// NOTE: yeah the confDB has the minimized energies already,
			// but it doesn't have the structures so we need to minimize again
			return analyzeEnsemble(table.energiedConfs(ConfDB.SortOrder.Energy).iterator(), maxNumConfs);
		});
	}

	public EnsembleAnalysis analyzeEnsemble(Queue.FIFO<? extends ConfSearch.ScoredConf> confs, int maxNumConfs) {
		return analyzeEnsemble(confs.iterator(), maxNumConfs);
	}

	public EnsembleAnalysis analyzeEnsemble(Iterator<? extends ConfSearch.ScoredConf> confs, int maxNumConfs) {

		EnsembleAnalysis analysis = new EnsembleAnalysis();

		// read the top confs
		for (int i=0; i<maxNumConfs; i++) {
			final int fi = i;

			// get the next conf
			if (!confs.hasNext()) {
				break;
			}
			ConfSearch.ScoredConf conf = confs.next();

			// minimize the conf (asynchronously if possible)
			confEcalc.tasks.submit(
				() -> confEcalc.calcEnergy(new RCTuple(conf.getAssignments())),
				(epmol) -> {
					while (analysis.analyses.size() <= fi) {
						analysis.analyses.add(null);
					}
					analysis.analyses.set(fi, new ConfAnalysis(
						conf.getAssignments(),
						conf.getScore(),
						epmol
					));
				}
			);
		}

		confEcalc.tasks.waitForFinish();

		return analysis;
	}
}
