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

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Cluster;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.AutoCloseableNoEx;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;


/**
 * Implementation of the K* algorithm to predict protein sequence mutations that improve
 * binding affinity by computing provably accurate Boltzmann-weighted ensembles
 * {@cite Lilien2008 Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson, and Bruce R. Donald, 2005.
 * A Novel Ensemble-Based Scoring and Search Algorithm for Protein Redesign and Its Application
 * to Modify the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine Adenylation Enzyme
 * In Journal of Computational Biology (vol 12. num. 6 pp. 740–761).}.
 */
public class KStar {

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		private final int maxNumConfs;

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
			 * True to use external memory when buffering conformations between the
			 * partition function lower and upper bound calculators.
			 */
			private boolean useExternalMemory = false;

			/**
			 * Pattern for ConfDB filename, where the first %s is replaced with the state name.
			 */
			private String confDBPattern = "%s.confdb";

			/**
			 * True to attempt to resume a previous design using the conformation databases.
			 * False to delete any existing conformation databases and start the design from scratch.
			 */
			private boolean resume = false;

			/**
			 * The maximum number of conformations for each pfunc to explore. Breaks provability, but
			 * can be used as a faster heuristic. -1 means use epsilon.
			 */
			private int maxNumberConfs = -1;

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

			public Builder setExternalMemory(boolean val) {
				useExternalMemory = val;
				return this;
			}

			public Builder setConfDBPattern(String val) {
				confDBPattern = val;
				return this;
			}

			public Builder setMaxNumConf(int val) {
				this.maxNumberConfs = val;
				return this;
			}

			public Builder resume(boolean val) {
				resume = val;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, stabilityThreshold, maxSimultaneousMutations, scoreWriters, showPfuncProgress, useExternalMemory, confDBPattern, resume, maxNumberConfs);
			}
		}

		public final double epsilon;
		public final Double stabilityThreshold;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;
		public final boolean showPfuncProgress;
		public final boolean useExternalMemory;
		public final String confDBPattern;
		public final boolean resume;

		public Settings(double epsilon, Double stabilityThreshold, int maxSimultaneousMutations, KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, boolean useExternalMemory, String confDBPattern, boolean resume, int maxNumberConfs) {
			this.epsilon = epsilon;
			this.stabilityThreshold = stabilityThreshold;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
			this.showPfuncProgress = dumpPfuncConfs;
			this.useExternalMemory = useExternalMemory;
			this.confDBPattern = confDBPattern;
			this.resume = resume;
			this.maxNumConfs = maxNumberConfs;
		}
	}

	public static class ScoredSequence {

		public final Sequence sequence;
		public final KStarScore score;

		public ScoredSequence(Sequence sequence, KStarScore score) {
			this.sequence = sequence;
			this.score = score;
		}

		@Override
		public String toString() {
			return "sequence: " + sequence + "   K*(log10): " + score;
		}

		public String toString(Sequence wildtype) {
			return "sequence: " + sequence.toString(Sequence.Renderer.AssignmentMutations) + "   K*(log10): " + score;
		}
	}

	public enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public static class InitException extends RuntimeException {

		public InitException(ConfSpaceType type, String name) {
			super(String.format("set %s for the %s conf space info before running", name, type.name()));
		}
	}

	public interface PfuncFactory {
		PartitionFunction make(RCs rcs);
	}

	public class ConfSpaceInfo {

		public final ConfSpaceIteration confSpace;
		public final ConfSpaceType type;
		public final String id;

		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public File confDBFile = null;
		public ConfEnergyCalculator confEcalc = null;
		public PfuncFactory pfuncFactory = null;

		private ConfDB confDB = null;

		public ConfSpaceInfo(ConfSpaceIteration confSpace, ConfSpaceType type) {
			this.confSpace = confSpace;
			this.type = type;
			this.id = type.name().toLowerCase();

			confDBFile = new File(String.format(settings.confDBPattern, id));
		}

		private void check() {
			if (confEcalc == null) {
				throw new InitException(type, "confEcalc");
			}
			if (pfuncFactory == null) {
				throw new InitException(type, "pfuncFactory");
			}
		}

		private AutoCloseableNoEx openConfDB() {
			if (confDBFile != null) {
				if (!settings.resume) {
					confDBFile.delete();
				}
				confDB = new ConfDB(confSpace, confDBFile);
			}
			return () -> {
				if (confDB != null) {
					confDB.close();
					confDB = null;
				}
			};
		}

		private PartitionFunction makePfunc(TaskExecutor.ContextGroup ctxGroup, Sequence seq) {

			RCs rcs = seq.makeRCs(confSpace);

			PartitionFunction pfunc = pfuncFactory.make(rcs);

			pfunc.setReportProgress(settings.showPfuncProgress);
			if (settings.useExternalMemory) {
				PartitionFunction.WithExternalMemory.setOrThrow(pfunc, true, rcs);
			}
			if (confDB != null) {
				PartitionFunction.WithConfDB.cast(pfunc).setConfDB(confDB, seq);
			}

			pfunc.setInstanceId(type.ordinal());

			pfunc.init(settings.epsilon);
			pfunc.putTaskContexts(ctxGroup);

			return pfunc;
		}

		public void clear() {
			pfuncResults.clear();
		}

		public PartitionFunction.Result calcPfunc(TaskExecutor.ContextGroup ctxGroup, Sequence globalSequence, BigDecimal stabilityThreshold) {

			Sequence sequence = globalSequence.filter(confSpace.seqSpace());

			// check the cache first
			PartitionFunction.Result result = pfuncResults.get(sequence);
			if (result != null) {
				return result;
			}

			// cache miss, need to compute the partition function

			// compute the partition function
			PartitionFunction pfunc = makePfunc(ctxGroup, sequence);
			pfunc.setStabilityThreshold(stabilityThreshold);
			pfunc.compute(settings.maxNumConfs);

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);

			/* HACKHACK: we're done using the A* tree, pfunc, etc
				and normally the garbage collector will clean them up,
				along with their off-heap resources (e.g. TPIE data structures).
				Except the garbage collector might not do it right away.
				If we try to allocate more off-heap resources before these get cleaned up,
				we might run out. So poke the garbage collector now and try to get
				it to clean up the off-heap resources right away.
			*/
			Runtime.getRuntime().gc();

			// newer JVMs have more concurrent garbage collectors
			// give it a little time to finish cleaning up the pfunc
			try {
				Thread.sleep(10);
			} catch (InterruptedException ex) {
				throw new RuntimeException(ex);
			}

			return result;
		}
	}

	private interface Scorer {
		KStarScore score(int sequenceNumber, PartitionFunction.Result proteinResult, PartitionFunction.Result ligandResult, PartitionFunction.Result complexResult);
	}

	/** A configuration space containing just the protein strand */
	public final ConfSpaceInfo protein;

	/** A configuration space containing just the ligand strand */
	public final ConfSpaceInfo ligand;

	/** A configuration space containing both the protein and ligand strands */
	public final ConfSpaceInfo complex;

	/** Optional and overridable settings for K* */
	public final Settings settings;

	private List<Sequence> sequences;

	public KStar(ConfSpaceIteration protein, ConfSpaceIteration ligand, ConfSpaceIteration complex, Settings settings) {
		this.settings = settings;
		this.protein = new ConfSpaceInfo(protein, ConfSpaceType.Protein);
		this.ligand = new ConfSpaceInfo(ligand, ConfSpaceType.Ligand);
		this.complex = new ConfSpaceInfo(complex, ConfSpaceType.Complex);
		this.sequences = new ArrayList<>();
	}

	public Iterable<ConfSpaceInfo> confSpaceInfos() {
		return Arrays.asList(protein, ligand, complex);
	}

	public ConfSpaceInfo getConfSpaceInfo(ConfSpaceIteration confSpace) {
		if (confSpace == protein.confSpace) {
			return protein;
		} else if (confSpace == ligand.confSpace) {
			return ligand;
		} else if (confSpace == complex.confSpace) {
			return complex;
		} else {
			throw new IllegalArgumentException("conf space does not match any known by this K* instance");
		}
	}

	public ScoredSequence score(Sequence seq, TaskExecutor tasks) {

		// make a context group for the task executor
		try (TaskExecutor.ContextGroup ctxGroup = tasks.contextGroup()) {

			// open the conf databases if needed
			try (AutoCloseableNoEx proteinCloser = protein.openConfDB()) {
			try (AutoCloseableNoEx ligandCloser = ligand.openConfDB()) {
			try (AutoCloseableNoEx complexCloser = complex.openConfDB()) {

				// check the conf space infos to make sure we have all the inputs
				protein.check();
				ligand.check();
				complex.check();

				// reset any previous state
				sequences.clear();
				protein.clear();
				ligand.clear();
				complex.clear();

				// skip the calculation on member nodes
				if (tasks instanceof Cluster.Member) {
					// TODO: try to get the scored sequence from the client?
					return null;
				}

				return new ScoredSequence(seq, new KStarScore(
					protein.calcPfunc(ctxGroup, seq, BigDecimal.ZERO),
					ligand.calcPfunc(ctxGroup, seq, BigDecimal.ZERO),
					complex.calcPfunc(ctxGroup, seq, BigDecimal.ZERO)
				));
			}}}
		}
	}

	public List<ScoredSequence> run() {
		// run without task contexts
		// useful for LUTE ecalcs, which don't use parallelism at all
		return run(new TaskExecutor());
	}

	public List<ScoredSequence> run(TaskExecutor tasks) {

		// make a context group for the task executor
		try (TaskExecutor.ContextGroup ctxGroup = tasks.contextGroup()) {

				// skip the calculation on member nodes
				if (tasks instanceof Cluster.Member) {
					// TODO: try to get the scored sequences from the client?
					return null;
				}

				// check the conf space infos to make sure we have all the inputs
				protein.check();
				ligand.check();
				complex.check();

				// reset any previous state
				sequences.clear();
				protein.clear();
				ligand.clear();
				complex.clear();


				List<ScoredSequence> scores = new ArrayList<>();

				// collect all the sequences explicitly
				if (complex.confSpace.seqSpace().containsWildTypeSequence()) {
					sequences.add(complex.confSpace.seqSpace().makeWildTypeSequence());
				}
				sequences.addAll(complex.confSpace.seqSpace().getMutants(settings.maxSimultaneousMutations, true));

				// TODO: sequence filtering? do we need to reject some mutation combinations for some reason?

				// now we know how many sequences there are in total
				int n = sequences.size();

				// make the sequence scorer and reporter
				Scorer scorer = (sequenceNumber, proteinResult, ligandResult, complexResult) -> {

					// compute the K* score
					KStarScore kstarScore = new KStarScore(proteinResult, ligandResult, complexResult);
					Sequence sequence = sequences.get(sequenceNumber);
					scores.add(new ScoredSequence(sequence, kstarScore));

					// report scores
					settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
						sequenceNumber,
						n,
						sequence,
						kstarScore,
						this
					));

					return kstarScore;
				};

				System.out.println("computing K* scores for " + sequences.size() + " sequences to epsilon = " + settings.epsilon + " ...");
				settings.scoreWriters.writeHeader();
				// TODO: progress bar?

				// open the conf databases if needed
				BigDecimal proteinStabilityThreshold = null;
				BigDecimal ligandStabilityThreshold = null;
				PartitionFunction.Result proteinResult;
				PartitionFunction.Result complexResult;
				PartitionFunction.Result ligandResult;

				try (AutoCloseableNoEx proteinCloser = protein.openConfDB()) {
				try (AutoCloseableNoEx ligandCloser = ligand.openConfDB()) {
				try (AutoCloseableNoEx complexCloser = complex.openConfDB()) {
					// compute wild type partition functions first (always at pos 0)
					proteinResult = protein.calcPfunc(ctxGroup, sequences.get(0), BigDecimal.ZERO);
					ligandResult = ligand.calcPfunc(ctxGroup, sequences.get(0), BigDecimal.ZERO);
					complexResult = complex.calcPfunc(ctxGroup, sequences.get(0), BigDecimal.ZERO);
				}}}

				KStarScore wildTypeScore = scorer.score(
						0,
						proteinResult,
						ligandResult,
						complexResult
				);

				if (settings.stabilityThreshold != null) {
					BigDecimal stabilityThresholdFactor = new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(settings.stabilityThreshold);
					proteinStabilityThreshold = wildTypeScore.protein.values.calcLowerBound().multiply(stabilityThresholdFactor);
					ligandStabilityThreshold = wildTypeScore.ligand.values.calcLowerBound().multiply(stabilityThresholdFactor);
				}

				// compute all the partition functions and K* scores for the rest of the sequences
				for (int i=1; i<n; i++) {

					try (AutoCloseableNoEx proteinCloser = protein.openConfDB()) {
					try (AutoCloseableNoEx ligandCloser = ligand.openConfDB()) {
					try (AutoCloseableNoEx complexCloser = complex.openConfDB()) {
						Sequence seq = sequences.get(i);

						// get the pfuncs, with short circuits as needed
						proteinResult = protein.calcPfunc(ctxGroup, seq, proteinStabilityThreshold);
						if (!KStarScore.isLigandComplexUseful(proteinResult)) {
							ligandResult = PartitionFunction.Result.makeAborted();
							complexResult = PartitionFunction.Result.makeAborted();
						} else {
							ligandResult = ligand.calcPfunc(ctxGroup, seq, ligandStabilityThreshold);
							if (!KStarScore.isComplexUseful(proteinResult, ligandResult)) {
								complexResult = PartitionFunction.Result.makeAborted();
							} else {
								complexResult = complex.calcPfunc(ctxGroup, seq, BigDecimal.ZERO);
							}
						}
					}}}

					scorer.score(i, proteinResult, ligandResult, complexResult);
				}

			return scores;
		}
	}
}
