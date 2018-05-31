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

package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.ThreadPoolTaskExecutor;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.*;
import smile.data.SparseDataset;
import smile.math.matrix.Matrix;
import smile.regression.LASSO;

import java.io.File;
import java.math.BigInteger;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;
import static edu.duke.cs.osprey.tools.Log.logf;


public class LUTE {

	public static enum Fitter {

		/**
		 * ordinary least squares via conjugate gradient
		 *
		 * ie, plain old linear regression
		 *
		 * TL;DR - OLSCG is faster than LASSO here, use OLSCG by default
		 *
		 * fastest to compute, tends to overfit a bit more than LASSO though, so requires more samples
		 * (even with the extra sample energies, CG is still generally faster overall than LASSO though)
		 */
		OLSCG(true) {

			@Override
			public double[] fit(LinearSystem system, LinearSystem.BInfo binfo, TaskExecutor tasks) {

				// build the linear model: Ax=b
				// except conjugate gradient needs square A, so transform to A^tAx = A^tb
				RealLinearOperator AtA = new RealLinearOperator() {

					@Override
					public int getRowDimension() {
						return system.tuples.size();
					}

					@Override
					public int getColumnDimension() {
						return system.tuples.size();
					}

					@Override
					public RealVector operate(RealVector vx)
							throws DimensionMismatchException {
						double[] x = ((ArrayRealVector)vx).getDataRef();
						double[] AtAx = system.parallelMultAt(system.parallelMultA(x, tasks), tasks);
						return new ArrayRealVector(AtAx, false);
					}
				};

				RealVector Atb = new ArrayRealVector(system.multAt(binfo.b), false);

				ConjugateGradient cg = new ConjugateGradient(100000, 1e-6, false);
				return ((ArrayRealVector)cg.solve(AtA, Atb)).getDataRef();
			}
		},

		/**
		 * least absolute shrinkage and selection operator
		 *
		 * ie, L1 regularized linear regression
		 *
		 * not as fast as conjugate gradient, but tends to overfit less,
		 * so usually requires fewer training samples. However, conjugate gradient
		 * is soooo much faster at fitting that it tends to be faster overall anyway.
		 */
		LASSO(false /* LASSO implementation does its own normalization */) {

			@Override
			public double[] fit(LinearSystem system, LinearSystem.BInfo binfo, TaskExecutor tasks) {

				// explicitly instantiate A (at least it's a sparse and not a dense matrix though)
				// the LASSO implementation is actually really fast!
				SparseDataset data = new SparseDataset(system.tuples.size());
				for (int c=0; c<system.confs.size(); c++) {
					final int fc = c;
					system.forEachTupleIn(c, (t) -> {
						data.set(fc, t, 1.0);
					});
				}
				Matrix A = data.toSparseMatrix();

				// tried a bunch of stuff, these seem to work well
				double lambda = 0.1;
				double tolerance = 0.1;
				int maxIterations = 200;

				// regress!
				// NOTE: can't do parallelism here apparently, so `tasks` is ingored =(
				LASSO lasso = new LASSO(A, binfo.b, lambda, tolerance, maxIterations);
				binfo.offset += lasso.intercept();
				return lasso.coefficients();
			}
		};

		public final boolean normalize;

		private Fitter(boolean normalize) {
			this.normalize = normalize;
		}

		public abstract double[] fit(LinearSystem system, LinearSystem.BInfo binfo, TaskExecutor tasks);
	}

	public static class LinearSystem {

		private class BInfo {
			double[] b = confEnergies.clone();
			double offset = 0.0;
			double scale = 1.0;
		}


		public final TuplesIndex tuples;
		public final List<int[]> confs;
		public final double[] confEnergies;

		public double[] tupleEnergies;
		public double tupleEnergyOffset;

		public Errors errors = null;

		public LinearSystem(TuplesIndex tuples, ConfSampler.Samples samples, Map<int[],Double> confEnergies) {

			this.tuples = tuples;

			// linearize the conformations
			confs = new ArrayList<>(samples.getAllConfs());

			// linearize the conf energies
			this.confEnergies = new double[confs.size()];
			for (int c=0; c<confs.size(); c++) {
				this.confEnergies[c] = confEnergies.get(confs.get(c));
			}
		}

		private void forEachTupleIn(int c, Consumer<Integer> callback) {
			final boolean throwIfMissingSingle = false; // we're not fitting singles
			final boolean throwIfMissingPair = true; // we always fit to dense pairs, confs shouldn't be using pruned pairs
			tuples.forEachIn(confs.get(c), throwIfMissingSingle, throwIfMissingPair, callback);
		}

		public void fit(Fitter fitter, TaskExecutor tasks) {

			// calculate b, and normalize if needed
			BInfo binfo = new BInfo();
			if (fitter.normalize) {

				// get the range of b
				double min = Double.POSITIVE_INFINITY;
				double max = Double.NEGATIVE_INFINITY;
				for (int c=0; c<confs.size(); c++) {
					double e = this.confEnergies[c];
					min = Math.min(min, e);
					max = Math.max(max, e);
				}

				binfo.offset = min;
				binfo.scale = max - min;

				for (int c=0; c<confs.size(); c++) {
					binfo.b[c] = (binfo.b[c] - binfo.offset)/binfo.scale;
				}
			}

			double[] x = fitter.fit(this, binfo, tasks);

			calcTupleEnergies(x, binfo);
		}

		private double[] multA(double[] x) {
			double[] out = new double[confs.size()];
			for (int c=0; c<confs.size(); c++) {
				final int fc = c;
				forEachTupleIn(c, (t) -> {
					out[fc] += x[t];
				});
			}
			return out;
		}

		private double[] multAt(double[] x) {
			double[] out = new double[tuples.size()];
			for (int c=0; c<confs.size(); c++) {
				double xc = x[c];
				forEachTupleIn(c, (t) -> {
					out[t] += xc;
				});
			}
			return out;
		}

		private double[] parallelMultA(double[] x, TaskExecutor tasks) {

			int numConfs = confs.size();
			int numThreads = tasks.getParallelism();
			int partitionSize = numConfs/numThreads;

			double[] out = new double[numConfs];

			// partition the confs equally among the threads
			for (int i=0; i<numThreads; i++) {
				int startC = i*partitionSize;
				int stopC;
				if (i == numThreads - 1) {
					stopC = numConfs;
				} else {
					stopC = (i+1)*partitionSize;
				}
				tasks.submit(
						() -> {
							for (int c=startC; c<stopC; c++) {
								final int fc = c;
								forEachTupleIn(fc, (t) -> {
									out[fc] += x[t];
								});
							}
							return null;
						},
						(ignored) -> {}
				);
			}
			tasks.waitForFinish();

			return out;
		}

		private double[] parallelMultAt(double[] x, TaskExecutor tasks) {

			int numConfs = confs.size();
			int numThreads = tasks.getParallelism();
			int partitionSize = numConfs/numThreads;

			int numTuples = tuples.size();
			double[] out = new double[numTuples];

			// partition the confs equally among the threads
			for (int i=0; i<numThreads; i++) {
				int startC = i*partitionSize;
				int stopC;
				if (i == numThreads - 1) {
					stopC = numConfs;
				} else {
					stopC = (i+1)*partitionSize;
				}
				tasks.submit(
						() -> {
							double[] threadOut = new double[numTuples];
							for (int c=startC; c<stopC; c++) {
								double xc = x[c];
								forEachTupleIn(c, (t) -> {
									threadOut[t] += xc;
								});
							}
							return threadOut;
						},
						(threadOut) -> {
							// aggregate thread results to the real out array
							for (int t=0; t<numTuples; t++) {
								out[t] += threadOut[t];
							}
						}
				);
			}
			tasks.waitForFinish();

			return out;
		}

		private void calcTupleEnergies(double[] x, BInfo binfo) {

			// calculate the tuple energies in un-normalized space
			double[] energies = new double[tuples.size()];
			for (int t=0; t<tuples.size(); t++) {
				energies[t] = x[t]*binfo.scale;
			}

			setTupleEnergies(energies, binfo.offset);
		}

		public void setTupleEnergies(double[] energies, double offset) {

			this.tupleEnergies = energies;
			this.tupleEnergyOffset = offset;

			// calculate the residual
			double[] residual = multA(tupleEnergies);
			for (int c=0; c<confs.size(); c++) {
				residual[c] = (residual[c] + tupleEnergyOffset - confEnergies[c]);
			}

			// analyze the errors
			this.errors = new Errors(residual);
		}
	}

	public static class Errors {

		public final double[] residual;
		public final double min;
		public final double max;
		public final double avg;
		public final double rms;

		public Errors(double[] residual) {

			this.residual = residual;

			double sum = 0.0;
			double sumsq = 0.0;
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;

			int n = residual.length;
			for (int row=0; row<n; row++) {

				double val = Math.abs(residual[row]);

				sum += val;
				sumsq += val*val;
				min = Math.min(min, val);
				max = Math.max(max, val);
			}

			this.min = min;
			this.max = max;
			this.avg = sum/n;
			this.rms = Math.sqrt(sumsq/n);
		}

		@Override
		public String toString() {
			return String.format("range [%8.4f,%8.4f]   avg %8.4f   rms %8.4f", min, max, avg, rms);
		}
	}

	public final SimpleConfSpace confSpace;

	private final Set<RCTuple> tuples = new LinkedHashSet<>();

	private LinearSystem trainingSystem = null;
	private LinearSystem testSystem = null;
	private Map<int[],Double> energies = null;

	public LUTE(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	public Set<RCTuple> getUnprunedSingleTuples(PruningMatrix pmat) {

		Set<RCTuple> singles = new LinkedHashSet<>();

		for (int pos1=0; pos1<pmat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<pmat.getNumConfAtPos(pos1); rc1++) {

				// skip pruned singles
				if (pmat.isSinglePruned(pos1, rc1)) {
					continue;
				}

				singles.add(new RCTuple(pos1, rc1));
			}
		}

		return singles;
	}

	public Set<RCTuple> getUnprunedPairTuples(PruningMatrix pmat) {

		Set<RCTuple> pairs = new LinkedHashSet<>();

		for (int pos1=0; pos1<pmat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<pmat.getNumConfAtPos(pos1); rc1++) {

				// skip pruned singles
				if (pmat.isSinglePruned(pos1, rc1)) {
					continue;
				}

				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<pmat.getNumConfAtPos(pos2); rc2++) {

						// skip pruned singles
						if (pmat.isSinglePruned(pos2, rc2)) {
							continue;
						}

						// skip pruned tuples
						if (pmat.isPairPruned(pos1, rc1, pos2, rc2)) {
							continue;
						}

						// we found it! It's an unpruned pair!
						// NOTE: make the tuple in pos2, pos1 order so the positions are already sorted
						// (because pos2 < pos1 by definition)
						pairs.add(new RCTuple(pos2, rc2, pos1, rc1));
					}
				}
			}
		}

		return pairs;
	}

	/**
	 * Turns out, this approach doesn't work terribly well unless you return lots of tuples
	 *
	 * I think it tends to overfit to the pairs training set, but doesn't fit well on later sampled confs for triples
	 */
	public Set<RCTuple> sampleTripleTuplesByFitError(PruningMatrix pmat, ConfDB.ConfTable confTable, double fractionSqErrorCovered) {

		LUTEConfEnergyCalculator luteConfEcalc = new LUTEConfEnergyCalculator(confSpace, new LUTEState(trainingSystem));

		class ConfInfo {

			public final int[] conf;
			public final double ffEnergy;
			public final double luteEnergy;
			public final double err;

			public ConfInfo(int[] conf) {
				this.conf = conf;
				ffEnergy = confTable.get(conf).upper.energy;
				luteEnergy = luteConfEcalc.calcEnergy(conf);
				err = Math.abs(ffEnergy - luteEnergy);
			}
		}

		// go through the conf table and try to find confs with big errors
		List<ConfInfo> confInfos = new ArrayList<>();
		for (int[] conf : trainingSystem.confs) {
			confInfos.add(new ConfInfo(conf));
		}
		confInfos.sort(Comparator.comparing((ConfInfo info) -> info.err).reversed());

		/* DEBUG
		log("%d training confs with energy gaps", confInfos.size());
		for (int i=0; i<100; i++) {

			ConfInfo info = confInfos.get(i);

			log("%30s:   ff %9.4f   lute %9.4f   err %9.4f",
				Conf.toString(info.conf),
				info.ffEnergy,
				info.luteEnergy,
				info.err
			);
		}
		*/

		// get all the triples we could possibly sample
		List<RCTuple> triples = new ArrayList<>();
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			triples.add(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted());
			return PruningMatrix.IteratorCommand.Continue;
		});


		// DEBUG
		//log("%d triples available", triples.size());

		// calculate an error score for each triple
		TupleTree<Double> triplesSqError = new TupleTree<>();
		for (ConfInfo info : confInfos) {

			for (RCTuple triple : triples) {
				if (Conf.containsTuple(info.conf, triple)) {

					Double errorSqSum = triplesSqError.get(triple);
					if (errorSqSum == null) {
						errorSqSum = 0.0;
					}
					errorSqSum += info.err*info.err;
					triplesSqError.put(triple, errorSqSum);
				}
			}
		}

		Function<RCTuple,Double> getSqError = (triple) -> {
			Double sqError = triplesSqError.get(triple);
			if (sqError == null) {
				sqError = 0.0;
			}
			return sqError;
		};

		/* DEBUG: find the triples with the most squared error
		triples.sort(Comparator.comparing(getSqError).reversed());
		for (int i=0; i<10; i++) {
			RCTuple triple = triples.get(i);
			log("%9.4f   %s", triplesSqError.get(triple), triple);
		}
		*/

		// sum the total error
		double totalSqError = 0.0;
		for (RCTuple triple : triples) {
			totalSqError += getSqError.apply(triple);
		}

		// take the tuples that represent the top X fraction of the squared error
		Set<RCTuple> chosenTriples = new HashSet<>();
		double sqError = 0.0;
		for (RCTuple triple : triples) {

			chosenTriples.add(triple);

			sqError += getSqError.apply(triple);
			if (sqError/totalSqError >= fractionSqErrorCovered) {
				break;
			}
		}

		// DEBUG
		//log("%d triples chosen", chosenTriples.size());

		return chosenTriples;
	}

	public Set<RCTuple> sampleTripleTuplesByStrongInteractions(EnergyMatrix emat, PruningMatrix pmat, int maxNumPairsPerPosition) {

		int n = emat.getNumPos();
		int numPairsPerPosition = Math.min(n - 1, maxNumPairsPerPosition);

		// find strongest pairwise energies for each pair of positions
		PosMatrixGeneric<Double> strongestInteractions = new PosMatrixGeneric<>(confSpace);
		strongestInteractions.fill(0.0);
		for (int pos1=1; pos1<n; pos1++) {
			for (int pos2=0; pos2<pos1; pos2++) {

				// get the strongest (ie, in magnitude only) energy over all RCs for this pair
				for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {

						double interaction = Math.abs(emat.getPairwise(pos1, rc1, pos2, rc2));

						if (interaction > strongestInteractions.get(pos1, pos2)) {
							strongestInteractions.set(pos1, pos2, interaction);
						}
					}
				}
			}
		}

		// pick the top k strongly-interacting positions for each position
		PosMatrixGeneric<Boolean> topPositionInteractions = new PosMatrixGeneric<>(confSpace);
		topPositionInteractions.fill(false);
		for (int pos1=0; pos1<n; pos1++) {

			// collect all positions that aren't pos1
			List<Integer> positions = new ArrayList<>(n);
			for (int pos2=0; pos2<n; pos2++) {
				if (pos2 != pos1) {
					positions.add(pos2);
				}
			}

			// sort them by strongest interactions with pos1
			final int fpos1 = pos1;
			positions.sort(Comparator.comparing((Integer pos2) -> topPositionInteractions.get(fpos1, pos2)).reversed());

			// pick the top k
			for (int i=0; i<numPairsPerPosition; i++) {
				int pos2 = positions.get(i);
				topPositionInteractions.set(pos1, pos2, true);
			}
		}

		// pick triples involving at least 2/3 strongly interacting pairs
		Set<RCTuple> triples = new HashSet<>();
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {

			int numStrongInteractions =
					(topPositionInteractions.get(pos1, pos2) ? 1 : 0)
							+ (topPositionInteractions.get(pos1, pos3) ? 1 : 0)
							+ (topPositionInteractions.get(pos2, pos3) ? 1 : 0);
			if (numStrongInteractions >= 2) {
				triples.add(new RCTuple(pos1, rc1, pos2, rc2, pos3, rc3).sorted());
			}

			return PruningMatrix.IteratorCommand.Continue;
		});
		return triples;
	}

	public void addTuples(Iterable<RCTuple> tuples) {
		for (RCTuple tuple : tuples) {
			addTuple(tuple);
		}
	}

	public void addTuple(RCTuple tuple) {
		tuple.checkSortedPositions();
		tuples.add(tuple);
	}

	public Set<RCTuple> getTuples() {
		// don't allow callers to directly modify the tuple set
		// we need to maintain the ordering of tuple positions
		return Collections.unmodifiableSet(tuples);
	}

	public void sampleTuplesAndFit(ConfEnergyCalculator confEcalc, EnergyMatrix emat, PruningMatrix pmat, ConfDB.ConfTable confTable, ConfSampler sampler, Fitter fitter, double maxOverfittingScore, double maxRMSE) {

		// does the conf space only have one position?
		if (confSpace.positions.size() == 1) {

			// can only do singles with such a small conf space
			logf("Sampling all pair tuples...");
			Stopwatch singlesStopwatch = new Stopwatch().start();
			addTuples(getUnprunedSingleTuples(pmat));
			log(" done in " + singlesStopwatch.stop().getTime(2));
			fit(confEcalc, confTable, sampler, fitter, maxOverfittingScore);

			// was that good enough?
			if (trainingSystem.errors.rms <= maxRMSE) {
				log("training set RMS error %f meets goal of %f", trainingSystem.errors.rms, maxRMSE);
			} else {
				log("training set RMS error %f does not meet goal of %f", trainingSystem.errors.rms, maxRMSE);
			}
			return;
		}

		// start with all pairs first
		logf("Sampling all pair tuples...");
		Stopwatch pairsStopwatch = new Stopwatch().start();
		addTuples(getUnprunedPairTuples(pmat));
		log(" done in " + pairsStopwatch.stop().getTime(2));
		fit(confEcalc, confTable, sampler, fitter, maxOverfittingScore);

		// was that good enough?
		if (trainingSystem.errors.rms <= maxRMSE) {
			log("training set RMS error %f meets goal of %f", trainingSystem.errors.rms, maxRMSE);
			return;
		}
		log("training set RMS error %f does not meet goal of %f", trainingSystem.errors.rms, maxRMSE);

		int maxNumPairsPerPosition = confSpace.positions.size() - 1;
		for (int numPairsPerPosition=1; numPairsPerPosition<maxNumPairsPerPosition; numPairsPerPosition++) {

			// nope, try to pick some triples to get a better fit
			logf("Sampling triple tuples (top %d strongly-interacting pairs at each position) to try to reduce error...", numPairsPerPosition);
			Stopwatch triplesStopwatch = new Stopwatch().start();
			addTuples(sampleTripleTuplesByStrongInteractions(emat, pmat, numPairsPerPosition));
			log(" done in " + triplesStopwatch.stop().getTime(2));

			// fit again
			fit(confEcalc, confTable, sampler, fitter, maxOverfittingScore);

			// was that good enough?
			if (trainingSystem.errors.rms <= maxRMSE) {
				log("training set RMS error %f meets goal of %f", trainingSystem.errors.rms, maxRMSE);
				return;
			}
			log("training set RMS error %f does not meet goal of %f", trainingSystem.errors.rms, maxRMSE);
		}

		log("all triples exhausted. Nothing more to try to improve the fit.");
	}

	public void fit(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, ConfSampler sampler, Fitter fitter, double maxOverfittingScore) {

		TuplesIndex tuplesIndex = new TuplesIndex(confSpace, tuples);

		ConfSampler.Samples trainingSet = new ConfSampler.Samples(tuplesIndex);
		ConfSampler.Samples testSet = new ConfSampler.Samples(tuplesIndex);
		energies = new Conf.Map<>();
		double overfittingScore;

		Stopwatch sw = new Stopwatch().start();

		int numSamples = 0;
		int samplesPerTuple = 1;
		while (true) {

			// sample the training and test sets
			logf("\nsampling at least %d confs per tuple for %d tuples...", samplesPerTuple, tuples.size());
			Stopwatch samplingSw = new Stopwatch().start();
			sampler.sampleConfsForTuples(trainingSet, samplesPerTuple);
			sampler.sampleConfsForTuples(testSet, samplesPerTuple);
			log(" done in %s", samplingSw.stop().getTime(2));

			// count all the samples
			for (int[] conf : trainingSet.getAllConfs()) {
				energies.put(conf, null);
			}
			for (int[] conf : testSet.getAllConfs()) {
				energies.put(conf, null);
			}

			int numAdditionalSamples = energies.size() - numSamples;
			numSamples = energies.size();

			// calculate energies for all the samples
			Progress progress = new Progress(numSamples);
			log("calculating energies for %d more samples...", numAdditionalSamples);
			for (Map.Entry<int[],Double> entry : energies.entrySet()) {
				int[] conf = entry.getKey();
				confEcalc.calcEnergyAsync(new RCTuple(conf), confTable, (energy) -> {
					entry.setValue(energy);
					progress.incrementProgress();
				});
			}
			confEcalc.tasks.waitForFinish();

			// make sure we get a thread pool intended for the CPU rather than the GPU
			try (ThreadPoolTaskExecutor tasks = new ThreadPoolTaskExecutor()) {
				tasks.start(confEcalc.ecalc.parallelism.numThreads);

				// fit the linear system to the training set
				logf("fitting %d confs to %d tuples ...", numSamples, this.tuples.size());
				Stopwatch trainingSw = new Stopwatch().start();
				trainingSystem = new LinearSystem(tuplesIndex, trainingSet, energies);
				trainingSystem.fit(fitter, tasks);
				logf(" done in %s", trainingSw.stop().getTime(2));
			}

			// analyze the test set errors
			testSystem = new LinearSystem(tuplesIndex, testSet, energies);
			testSystem.setTupleEnergies(trainingSystem.tupleEnergies, trainingSystem.tupleEnergyOffset);

			// measure overfitting by comparing ratio of rms errors
			overfittingScore = calcOverfittingScore();
			if (overfittingScore <= maxOverfittingScore) {
				log("");
				break;
			}

			log("    RMS errors:  train %.4f    test %.4f    overfitting score: %.4f",
					trainingSystem.errors.rms,
					testSystem.errors.rms,
					overfittingScore
			);

			samplesPerTuple++;
		}

		log("\nLUTE fitting finished in %s:\n", sw.stop().getTime(2));
		log("sampled at least %d confs per tuple for %d tuples", samplesPerTuple, tuples.size());
		log("sampled %d training confs, %d test confs, %d confs total (%.2f%% overlap)",
				trainingSet.size(), testSet.size(), energies.size(),
				100.0*(trainingSet.size() + testSet.size() - energies.size())/energies.size()
		);
		log("training errors: %s", trainingSystem.errors);
		log("    test errors: %s", testSystem.errors);

		// compute a histogram of test set errors
		double[] bucketTops = { 0.001, 0.01, 0.1, 0.2, 0.4, 1.0, 2.0, 4.0, 10.0 };
		int[] counts = new int[bucketTops.length];
		Arrays.fill(counts, 0);
		for (double error : testSystem.errors.residual) {
			for (int i=0; i<counts.length; i++) {
				if (Math.abs(error) <= bucketTops[i]) {
					counts[i]++;
				}
			}
		}
		for (int i=0; i<counts.length; i++) {
			log("               : %5.1f%% <= %s", 100.0*counts[i]/testSystem.confs.size(), Double.toString(bucketTops[i]));
		}

		log("total energy calculations: %d    overfitting score: %.4f <= %.4f", energies.size(), overfittingScore, maxOverfittingScore);
		log("");
	}

	public double calcOverfittingScore() {

		double num = testSystem.errors.rms;
		double denom = trainingSystem.errors.rms;

		if (num == 0.0 && denom == 0.0) {
			// technically undefined, but let's assume we're not overfitting at all in this case
			return 1.0;
		} else {
			return num/denom;
		}
	}

	public void reportConfSpaceSize() {

		BigInteger size = new RCs(confSpace).getNumConformations();

		double percent = 100.0*energies.size()/size.doubleValue();

		log("conf space (no pruning was reported) has exactly %s conformations", formatBig(size));
		log("LUTE sampled %.1f percent of those confs", percent);
	}

	public void reportConfSpaceSize(PruningMatrix pmat) {

		BigInteger sizeLower = pmat.calcUnprunedConfsLowerBound();
		BigInteger sizeUpper = pmat.calcUnprunedConfsUpperBound();

		// the bounds are loose, but remember we sampled some confs
		// so improve the lower bound with the samples if possible
		try {
			sizeLower = BigInteger.valueOf(Math.max(sizeLower.longValueExact(), energies.size()));

		} catch (ArithmeticException ex) {
			// lower bound is bigger than we could have possibly sampled
			// so don't change anything
		}

		double percentLower = 100.0*energies.size()/sizeUpper.doubleValue();
		double percentUpper = 100.0*energies.size()/sizeLower.doubleValue();

		log("conf space (after singles and pairs pruning) has somewhere between %s and %s conformations", formatBig(sizeLower), formatBig(sizeUpper));
		log("LUTE sampled somewhere between %.1f%% and %.1f%% of those conformations", percentLower, percentUpper);
	}

	public LinearSystem getTrainingSystem() {
		return trainingSystem;
	}

	public LinearSystem getTestSystem() {
		return testSystem;
	}

	public void save(File file) {
		LUTEIO.write(new LUTEState(trainingSystem), file);
	}
}
