package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Implementation of the K* algorithm to predict protein sequence mutations that improve
 * binding affinity by computing provably accurate Boltzmann-weighted ensembles
 * {@cite Lilien2008 Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson, and Bruce R. Donald, 2005.
 * A Novel Ensemble-Based Scoring and Search Algorithm for Protein Redesign and Its Application
 * to Modify the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine Adenylation Enzyme
 * In Journal of Computational Biology (vol 12. num. 6 pp. 740–761).}.
 */
public class KStar {

	public static interface ConfEnergyCalculatorFactory {
		ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	public static interface ConfSearchFactory {
		public ConfSearch make(EnergyMatrix emat, RCs rcs);
	}

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			/**
			 * Value of epsilon in (0,1] for the epsilon-approximation to a partition function.
			 *
			 * Smaller values for epsilon yield more accurate predictions, but can take
			 * longer to run.
			 */
			private double epsilon = 0.683;

			/**
			 * Pruning criteria to remove sequences with unstable unbound states relative to the wild type sequence.
			 * Defined in units of kcal/mol.
			 *
			 * More precisely, a sequence is pruned when the following expression is true:
			 *
			 * U(Z_s) < L(W_s) * B(t)
			 *
			 * where:
			 *   - s represents the unbound protein strand, or unbound ligand strand
			 *   - U(Z_s) is the upper bound on the partition function for strand s
			 *   - L(W_s) is the lower bound on the partition function for strand s in the wild type
			 *   - t is the stability threshold
			 *   - B() is the Boltzmann weighting function
			 *
			 * Set to null to disable the filter entirely.
			 */
			private Double stabilityThreshold = 5.0;

			/** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
			private int maxSimultaneousMutations = 1;

			private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();

			/**
			 * If true, prints out information to the console for each minimized conformation during
			 * partition function approximation
			 */
			private boolean showPfuncProgress = false;

			/**
			 * Pattern of the filename to cache energy matrices.
			 *
			 * K*-type algorithms must calculate multiple energy matrices.
			 * By default, these energy matrices are not cached between runs.
			 * To cache energy matrices between runs, supply a pattern such as:
			 *
			 * "theFolder/emat.*.dat"
			 *
			 * The * in the pattern is a wildcard character that will be replaced with
			 * each type of energy matrix used by the K*-type algorithm.
			 */
			private String energyMatrixCachePattern = null;

			public Builder setEpsilon(double val) {
				epsilon = val;
				return this;
			}

			public Builder setStabilityThreshold(Double val) {
				if (val != null && val.isInfinite()) {
					throw new IllegalArgumentException("only finite values allowed. To turn off the filter, pass null");
				}
				stabilityThreshold = val;
				return this;
			}

			public Builder setMaxSimultaneousMutations(int val) {
				maxSimultaneousMutations = val;
				return this;
			}

			public Builder addScoreWriter(KStarScoreWriter val) {
				scoreWriters.add(val);
				return this;
			}

			public Builder addScoreConsoleWriter(KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToConsole(val));
			}

			public Builder addScoreConsoleWriter() {
				return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequenceKStarPfuncs());
			}

			public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
			}

			public Builder addScoreFileWriter(File file) {
				return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
			}

			public Builder setShowPfuncProgress(boolean val) {
				showPfuncProgress = val;
				return this;
			}

			public Builder setEnergyMatrixCachePattern(String val) {
				energyMatrixCachePattern = val;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, stabilityThreshold, maxSimultaneousMutations, scoreWriters, showPfuncProgress, energyMatrixCachePattern);
			}
		}

		public final double epsilon;
		public final Double stabilityThreshold;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;
		public final boolean showPfuncProgress;
		public final String energyMatrixCachePattern;

		public Settings(double epsilon, Double stabilityThreshold, int maxSimultaneousMutations, KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, String energyMatrixCachePattern) {
			this.epsilon = epsilon;
			this.stabilityThreshold = stabilityThreshold;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
			this.showPfuncProgress = dumpPfuncConfs;
			this.energyMatrixCachePattern = energyMatrixCachePattern;
		}

		public String applyEnergyMatrixCachePattern(String type) {

			// the pattern has a * right?
			if (energyMatrixCachePattern.indexOf('*') < 0) {
				throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + energyMatrixCachePattern + "') has no wildcard character (which is *)");
			}

			return energyMatrixCachePattern.replace("*", type);
		}
	}

	public static class Sequence extends ArrayList<String> {

		public static Sequence makeWildType(SimpleConfSpace confSpace) {
			return new KStar.Sequence(confSpace.positions.size()).fillWildType(confSpace);
		}

		public Sequence(int size) {
			super(size);
			for (int i=0; i<size; i++) {
				add(null);
			}
		}

		public Sequence(String resTypes) {
			this(Arrays.asList(resTypes.split(" ")));
		}

		public Sequence(List<String> resTypes) {
			super(resTypes);
		}

		public Sequence makeWithAssignment(int posIndex, String resType) {
			Sequence assigned = new Sequence(this);
			assigned.set(posIndex, resType);
			return assigned;
		}

		public Sequence fillWildType(SimpleConfSpace confSpace) {
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				if (get(pos.index) == null) {
					set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
				}
			}
			return this;
		}

		public int countAssignments() {
			int count = 0;
			for (String resType : this) {
				if (resType != null) {
					count++;
				}
			}
			return count;
		}

		public boolean isFullyAssigned() {
			for (String resType : this) {
				if (resType == null) {
					return false;
				}
			}
			return true;
		}

		public int countMutations(SimpleConfSpace confSpace) {
			int count = 0;
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				String resType = get(pos.index);
				if (resType != null && !resType.equals(pos.resFlex.wildType)) {
					count++;
				}
			}
			return count;
		}

		public boolean isWildType(SimpleConfSpace confSpace) {
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				String resType = get(pos.index);
				if (resType != null && !resType.equals(pos.resFlex.wildType)) {
					return false;
				}
			}
			return true;
		}

		public RCs makeRCs(SimpleConfSpace confSpace) {
			return new RCs(confSpace, (pos, resConf) -> {

				// if there's an assignment here, only keep matching RCs
				String resType = get(pos.index);
				if (resType != null) {
					return resConf.template.name.equalsIgnoreCase(resType);
				}

				// otherwise, keep everything
				return true;
			});
		}

		@Override
		public String toString() {
			return String.join(" ", this);
		}

		/** highlight mutations in upper case, show wild type residues in lower case */
		public String toString(KStar.Sequence wildtype) {
			return toString(wildtype, 3);
		}

		public String toString(KStar.Sequence wildtype, int cellSize) {
			List<String> resTypes = new ArrayList<>();
			for (int i=0; i<size(); i++) {
				String resType = get(i);
				String resWildType = wildtype.get(i);
				if (resType == null) {
					resTypes.add(null);
				} else if (resType.equalsIgnoreCase(resWildType)) {
					resTypes.add(resType.toLowerCase());
				} else {
					resTypes.add(resType.toUpperCase());
				}
			}
			return String.join(" ", resTypes.stream()
				.map((s) -> String.format("%-" + cellSize + "s", s))
				.collect(Collectors.toList())
			);
		}

		public String toString(List<SimpleConfSpace.Position> positions) {
			return String.join(" ", positions.stream()
				.map((pos) -> this.get(pos.index) + "-" + pos.resNum)
				.collect(Collectors.toList())
			);
		}
	}

	public static class ScoredSequence {

		public final KStar.Sequence sequence;
		public final KStarScore score;

		public ScoredSequence(KStar.Sequence sequence, KStarScore score) {
			this.sequence = sequence;
			this.score = score;
		}

		@Override
		public String toString() {
			return "sequence: " + sequence + "   K*(log10): " + score;
		}

		public String toString(KStar.Sequence wildtype) {
			return "sequence: " + sequence.toString(wildtype) + "   K*(log10): " + score;
		}
	}

	public static enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public class ConfSpaceInfo {

		public final ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;

		public final List<Sequence> sequences = new ArrayList<>();
		public EnergyMatrix emat = null;
		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc) {
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

		public Sequence makeWildTypeSequence() {
			Sequence sequence = new Sequence(confSpace.positions.size());
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				sequence.set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
			}
			return sequence;
		}

		public PartitionFunction.Result calcPfunc(int sequenceIndex, BigDecimal stabilityThreshold) {

			Sequence sequence = sequences.get(sequenceIndex);

			// check the cache first
			PartitionFunction.Result result = pfuncResults.get(sequence);
			if (result != null) {
				return result;
			}

			// cache miss, need to compute the partition function

			// make the partition function
			RCs rcs = sequence.makeRCs(confSpace);
			ConfSearch astar = confSearchFactory.make(emat, rcs);
			PartitionFunction pfunc = new SimplePartitionFunction(astar, confEcalc);
			pfunc.setReportProgress(settings.showPfuncProgress);

			// compute it
			pfunc.init(settings.epsilon, stabilityThreshold);
			pfunc.compute();

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);
			return result;
		}
	}

	private static interface Scorer {
		KStarScore score(int sequenceNumber, PartitionFunction.Result proteinResult, PartitionFunction.Result ligandResult, PartitionFunction.Result complexResult);
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
	public final ConfEnergyCalculatorFactory confEcalcFactory;

	/** A function that makes a ConfSearchFactory (e.g, A* search) with the desired options */
	public final ConfSearchFactory confSearchFactory;

	/** Optional and overridable settings for K* */
	public final Settings settings;

	public KStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, ConfEnergyCalculatorFactory confEcalcFactory, ConfSearchFactory confSearchFactory, Settings settings) {
		this.protein = new ConfSpaceInfo(ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc));
		this.ligand = new ConfSpaceInfo(ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc));
		this.complex = new ConfSpaceInfo(ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc));
		this.ecalc = ecalc;
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;
	}

	public List<ScoredSequence> run() {

		List<ScoredSequence> scores = new ArrayList<>();

		// compute energy matrices
		protein.calcEmat();
		ligand.calcEmat();
		complex.calcEmat();

		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToProteinMap = complex.confSpace.mapPositionsTo(protein.confSpace);
		Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToLigandMap = complex.confSpace.mapPositionsTo(ligand.confSpace);

		// collect the wild type sequences
		protein.sequences.add(protein.makeWildTypeSequence());
		ligand.sequences.add(ligand.makeWildTypeSequence());
		complex.sequences.add(complex.makeWildTypeSequence());

		// collect all the sequences explicitly
		List<List<SimpleConfSpace.Position>> powersetOfPositions = MathTools.powersetUpTo(complex.confSpace.positions, settings.maxSimultaneousMutations);
		Collections.reverse(powersetOfPositions); // NOTE: reverse to match order of old code
		for (List<SimpleConfSpace.Position> mutablePositions : powersetOfPositions) {

			// collect the mutations (res types except for wild type) for these positions into a simple list list
			List<List<String>> resTypes = new ArrayList<>();
			for (SimpleConfSpace.Position pos : mutablePositions) {
				resTypes.add(pos.resFlex.resTypes.stream()
					.filter((resType) -> !resType.equals(pos.resFlex.wildType))
					.collect(Collectors.toList())
				);
			}

			// enumerate all the combinations of res types
			for (List<String> mutations : MathTools.cartesianProduct(resTypes)) {

				// build the complex sequence
				Sequence complexSequence = complex.makeWildTypeSequence();
				for (int i=0; i<mutablePositions.size(); i++) {
					complexSequence.set(mutablePositions.get(i).index, mutations.get(i));
				}
				complex.sequences.add(complexSequence);

				// split complex sequence into protein/ligand sequences
				Sequence proteinSequence = protein.makeWildTypeSequence();
				Sequence ligandSequence = ligand.makeWildTypeSequence();
				for (SimpleConfSpace.Position pos : complex.confSpace.positions) {

					SimpleConfSpace.Position proteinPos = complexToProteinMap.get(pos);
					if (proteinPos != null) {
						proteinSequence.set(proteinPos.index, complexSequence.get(pos.index));
					}

					SimpleConfSpace.Position ligandPos = complexToLigandMap.get(pos);
					if (ligandPos != null) {
						ligandSequence.set(ligandPos.index, complexSequence.get(pos.index));
					}
				}
				protein.sequences.add(proteinSequence);
				ligand.sequences.add(ligandSequence);
			}
		}

		// TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

		// now we know how many sequences there are in total
		int n = complex.sequences.size();

		// make the sequence scorer and reporter
		Scorer scorer = (sequenceNumber, proteinResult, ligandResult, complexResult) -> {

			// compute the K* score
			KStarScore kstarScore = new KStarScore(proteinResult, ligandResult, complexResult);
			Sequence complexSequence = complex.sequences.get(sequenceNumber);
			scores.add(new ScoredSequence(complexSequence, kstarScore));

			// report scores
			settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
				sequenceNumber,
				n,
				complexSequence,
				complex.confSpace,
				kstarScore
			));

			return kstarScore;
		};

		System.out.println("computing K* scores for " + complex.sequences.size() + " sequences to epsilon = " + settings.epsilon + " ...");
		settings.scoreWriters.writeHeader();
		// TODO: progress bar?

		// compute wild type partition functions first (always at pos 0)
		KStarScore wildTypeScore = scorer.score(
			0,
			protein.calcPfunc(0, BigDecimal.ZERO),
			ligand.calcPfunc(0, BigDecimal.ZERO),
			complex.calcPfunc(0, BigDecimal.ZERO)
		);
		BigDecimal proteinStabilityThreshold = null;
		BigDecimal ligandStabilityThreshold = null;
		if (settings.stabilityThreshold != null) {
			BigDecimal stabilityThresholdFactor = new BoltzmannCalculator().calc(settings.stabilityThreshold);
			proteinStabilityThreshold = wildTypeScore.protein.values.calcLowerBound().multiply(stabilityThresholdFactor);
			ligandStabilityThreshold = wildTypeScore.ligand.values.calcLowerBound().multiply(stabilityThresholdFactor);
		}

		// compute all the partition functions and K* scores for the rest of the sequences
		for (int i=1; i<n; i++) {

			// get the pfuncs, with short circuits as needed
			final PartitionFunction.Result proteinResult = protein.calcPfunc(i, proteinStabilityThreshold);
			final PartitionFunction.Result ligandResult;
			final PartitionFunction.Result complexResult;
			if (!KStarScore.isLigandComplexUseful(proteinResult)) {
				ligandResult = PartitionFunction.Result.makeAborted();
				complexResult = PartitionFunction.Result.makeAborted();
			} else {
				ligandResult = ligand.calcPfunc(i, ligandStabilityThreshold);
				if (!KStarScore.isComplexUseful(proteinResult, ligandResult)) {
					complexResult = PartitionFunction.Result.makeAborted();
				} else {
					complexResult = complex.calcPfunc(i, BigDecimal.ZERO);
				}
			}

			scorer.score(i, proteinResult, ligandResult, complexResult);
		}

		return scores;
	}
}
