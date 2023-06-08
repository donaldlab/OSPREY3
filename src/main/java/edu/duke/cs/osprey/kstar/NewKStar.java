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

import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.AutoCloseableNoEx;

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
public class NewKStar {

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
	public final KStarSettings settings;

	private final List<Sequence> sequences;

	public NewKStar(ConfSpaceIteration protein, ConfSpaceIteration ligand, ConfSpaceIteration complex, KStarSettings settings) {
		this.settings = settings;
		this.protein = new ConfSpaceInfo(settings, protein, ConfSpaceType.Protein);
		this.ligand = new ConfSpaceInfo(settings, ligand, ConfSpaceType.Ligand);
		this.complex = new ConfSpaceInfo(settings, complex, ConfSpaceType.Complex);
		this.sequences = new ArrayList<>();
	}

	public List<ConfSpaceInfo> confSpaceInfos() {
		return Arrays.asList(protein, ligand, complex);
	}

	public ScoredSequence score(Sequence seq, TaskExecutor tasks) {

		// make a context group for the task executor
		try (var ctxGroup = tasks.contextGroup()) {
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
		try(var ctxGroup = tasks.contextGroup()) {

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

				System.out.printf("Computing sequence %d%n", i);
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
