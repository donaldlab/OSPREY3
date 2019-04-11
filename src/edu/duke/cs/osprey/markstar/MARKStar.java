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

package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStarScore;
import edu.duke.cs.osprey.kstar.KStarScoreWriter;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundFastQueues;
import edu.duke.cs.osprey.markstar.framework.MARKStarBoundRigid;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.Stopwatch;

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
public class MARKStar {

	public interface ConfEnergyCalculatorFactory {
		ConfEnergyCalculator make(SimpleConfSpace confSpace, EnergyCalculator ecalc);
	}

	public interface ConfSearchFactory {
		ConfSearch make(EnergyMatrix emat, RCs rcs);
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

			private Parallelism parallelism = null;
			private int maxNumConfs = -1;
			private boolean reduceMinimizations = true;

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

			public Builder setParallelism(Parallelism p) {
				parallelism = p;
				return this;
			}

			public Builder setMaxNumConfs(int maxNumConfs) {
				this.maxNumConfs = maxNumConfs;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, stabilityThreshold, maxSimultaneousMutations, scoreWriters,
						showPfuncProgress, energyMatrixCachePattern, parallelism, maxNumConfs, reduceMinimizations);
			}

			public Builder setReduceMinimizations(boolean reudceMinimizations) {
			    this.reduceMinimizations = reudceMinimizations;
			    return this;
			}
		}

		public final double epsilon;
		public final Double stabilityThreshold;
		public final int maxSimultaneousMutations;
		public final KStarScoreWriter.Writers scoreWriters;
		public final boolean showPfuncProgress;
		public final String energyMatrixCachePattern;
		public final Parallelism parallelism;
		public final int maxNumConfs;
		public final boolean reduceMinimizations;

		public Settings(double epsilon, Double stabilityThreshold, int maxSimultaneousMutations,
						KStarScoreWriter.Writers scoreWriters, boolean dumpPfuncConfs, String energyMatrixCachePattern,
						Parallelism parallelism, int maxNumConfs, boolean reduceMinimizations) {
			this.epsilon = epsilon;
			this.stabilityThreshold = stabilityThreshold;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.scoreWriters = scoreWriters;
			this.showPfuncProgress = dumpPfuncConfs;
			this.energyMatrixCachePattern = energyMatrixCachePattern;
			this.parallelism = parallelism;
			this.maxNumConfs = maxNumConfs;
			this.reduceMinimizations = reduceMinimizations;
		}

		public String applyEnergyMatrixCachePattern(String type) {

			// the pattern has a * right?
			if (energyMatrixCachePattern.indexOf('*') < 0) {
				throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + energyMatrixCachePattern + "') has no wildcard character (which is *)");
			}

			return energyMatrixCachePattern.replace("*", type);
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

	public static class InitException extends RuntimeException {

		public InitException(ConfSpaceType type, String name) {
			super(String.format("set %s for the %s conf space info before running", name, type.name()));
		}
	}
	public enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public class ConfSpaceInfo {

		public final ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator rigidConfEcalc;
		public final ConfEnergyCalculator minimizingConfEcalc;
		public UpdatingEnergyMatrix correctionEmat;
		public ConfSearchFactory confSearchFactory = null;
		public File confDBFile = null;

		public EnergyMatrix rigidEmat = null;
		public EnergyMatrix minimizingEmat = null;
		public final Map<Sequence,PartitionFunction.Result> pfuncResults = new HashMap<>();

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator rigidConfEcalc, ConfEnergyCalculator minimizingConfEcalc) {
			this.type = type;
			this.confSpace = confSpace;
			this.rigidConfEcalc = rigidConfEcalc;
			this.minimizingConfEcalc = minimizingConfEcalc;
		}

		private void check() {
			if (rigidConfEcalc == null) {
				throw new InitException(type, "rigidConfEcalc");
			}
			if (minimizingConfEcalc == null) {
				throw new InitException(type, "minimizingConfEcalc");
			}
			if (confSearchFactory == null) {
				throw new InitException(type, "confSearchFactory");
			}
		}

		public void clear() {
			pfuncResults.clear();
		}

		public void calcEmats() {
			SimplerEnergyMatrixCalculator.Builder rigidBuilder = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc);
			if (settings.energyMatrixCachePattern != null) {
				rigidBuilder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase()+".rigid")));
			}
			SimplerEnergyMatrixCalculator.Builder minimizingBuilder = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc);
			if (settings.energyMatrixCachePattern != null) {
				minimizingBuilder.setCacheFile(new File(settings.applyEnergyMatrixCachePattern(type.name().toLowerCase()+".minimizing")));
			}
			rigidEmat = rigidBuilder.build().calcEnergyMatrix();
			minimizingEmat = minimizingBuilder.build().calcEnergyMatrix();
			correctionEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat);
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
			MARKStarBoundFastQueues pfunc = new MARKStarBoundFastQueues(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, sequence.makeRCs(confSpace),
			//MARKStarBoundRigid pfunc = new MARKStarBoundRigid(confSpace, rigidEmat, minimizingEmat, minimizingConfEcalc, sequence.makeRCs(confSpace),
					settings.parallelism);
			confSearchFactory = (emat, rcs) -> {
				ConfAStarTree.Builder builder = new ConfAStarTree.Builder(emat, rcs)
						.setTraditional();
				return builder.build();
			};
			RCs rcs = sequence.makeRCs(confSpace);
			ConfSearch astar = confSearchFactory.make(minimizingEmat, rcs);
			//GradientDescentMARKStarPfunc pfunc = new GradientDescentMARKStarPfunc(confSpace, rigidEmat, minimizingEmat,
			//		rcs, minimizingConfEcalc);
			pfunc.reduceMinimizations = settings.reduceMinimizations;
			pfunc.stateName = type.name();
			if(settings.maxNumConfs > 0)
				pfunc.setMaxNumConfs(settings.maxNumConfs);
			pfunc.setReportProgress(settings.showPfuncProgress);

			pfunc.setCorrections(correctionEmat);

			if (settings.showPfuncProgress == true){
				System.out.println("Computing "+type+":");
			}

			// compute it
			pfunc.init(astar, rcs.getNumConformations(), settings.epsilon);
			Stopwatch computeTimer = new Stopwatch().start();
			pfunc.compute();
			computeTimer.stop();
			System.out.println("Computation for "+sequence.toString()+":"+computeTimer.getTime(2));

			// save the result
			result = pfunc.makeResult();
			pfuncResults.put(sequence, result);
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

	/** Calculates the rigid energy for a molecule */
	public final EnergyCalculator rigidEcalc;

	/** Calculates the minimized energy for a molecule */
	public final EnergyCalculator minimizingEcalc;

	/** A function that makes a ConfEnergyCalculator with the desired options */
	public final ConfEnergyCalculatorFactory confEcalcFactory;

	/** Optional and overridable settings for K* */
	public final Settings settings;

	private List<Sequence> sequences;

	public MARKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex,
					EnergyCalculator rigidEcalc, EnergyCalculator minimizingEcalc,
					ConfEnergyCalculatorFactory confEcalcFactory, Settings settings) {
		this.protein = new ConfSpaceInfo(ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, rigidEcalc), confEcalcFactory.make(protein, minimizingEcalc));
		this.ligand = new ConfSpaceInfo(ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, rigidEcalc), confEcalcFactory.make(ligand, minimizingEcalc));
		this.complex = new ConfSpaceInfo(ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, rigidEcalc), confEcalcFactory.make(complex, minimizingEcalc));
		this.rigidEcalc = rigidEcalc;
		this.minimizingEcalc = minimizingEcalc;
		this.confEcalcFactory = confEcalcFactory;
		this.settings = settings;
		this.sequences = new ArrayList();
	}

	public void precalcEmats() {
		// compute energy matrices
		protein.calcEmats();
		ligand.calcEmats();
		complex.calcEmats();
	}

	public List<ScoredSequence> run() {

		List<ScoredSequence> scores = new ArrayList<>();

		// compute energy matrices
		protein.calcEmats();
		ligand.calcEmats();
		complex.calcEmats();


		// collect all the seque// collect all the sequences explicitly
		if (complex.confSpace.seqSpace.containsWildTypeSequence()) {
			sequences.add(complex.confSpace.seqSpace.makeWildTypeSequence());
		}
		sequences.addAll(complex.confSpace.seqSpace.getMutants(settings.maxSimultaneousMutations, true));

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
				kstarScore
			));

			return kstarScore;
		};

		System.out.println("computing K* scores for " + sequences.size() + " sequences to epsilon = " + settings.epsilon + " ...");
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
			BigDecimal stabilityThresholdFactor = new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(settings.stabilityThreshold);
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
	public Iterable<ConfSpaceInfo> confSpaceInfos() {
		return Arrays.asList(protein, ligand, complex);
	}
}
