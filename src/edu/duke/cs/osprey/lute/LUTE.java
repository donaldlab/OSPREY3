package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfDB;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.Stopwatch;
import org.apache.commons.math3.linear.*;

import java.math.BigInteger;
import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;


public class LUTE {

	private static class LinearSystem {

		final List<RCTuple> tuples;
		final List<int[]> confs;
		final RealMatrix A;
		final RealVector b;

		public LinearSystem(Set<RCTuple> tuples, Map<RCTuple,Set<int[]>> samplesByTuple, Map<int[],Double> energies) {

			// linearize the tuples
			this.tuples = new ArrayList<>(tuples);

			// linearize the conformations
			Set<int[]> confsSet = new Conf.Set();
			for (Collection<int[]> samples : samplesByTuple.values()) {
				confsSet.addAll(samples);
			}
			confs = new ArrayList<>(confsSet);

			// build the linear model: Ax=b
			A = new OpenMapRealMatrix(confs.size(), this.tuples.size());
			b = new ArrayRealVector(confs.size());
			for (int row=0; row<confs.size(); row++) {
				int[] conf = confs.get(row);
				for (int col=0; col<this.tuples.size(); col++) {
					RCTuple tuple = this.tuples.get(col);
					if (Conf.containsTuple(conf, tuple)) {
						A.setEntry(row, col, 1);
					}
				}
				b.setEntry(row, energies.get(conf));
			}

			// TEMP
			log("confs: %d   tuples: %d", confs.size(), this.tuples.size());
			double minVal = Double.POSITIVE_INFINITY;
			double maxVal = Double.NEGATIVE_INFINITY;
			for (int row=0; row<confs.size(); row++) {
				double energy = b.getEntry(row);
				minVal = Math.min(minVal, energy);
				maxVal = Math.max(maxVal, energy);
			}
			log("energy range: [%.4f,%.4f]", minVal, maxVal);
		}

		public RealVector fit() {

			log("solving linear model...");

			// any linear least squares solver should work here
			// I generally like to use QR factorization, but it's a bit slow here
			// so let's use conjugate gradient iterative optimization instead

			//return new QRDecomposition(A).getSolver().solve(b);

			// need square A though, so transform Ax=b to A^tAx = A^tb
			RealMatrix At = A.transpose();
			RealMatrix AtA = At.multiply(A);
			RealVector Atb = At.operate(b);
			return new ConjugateGradient(100000, 1e-6, false).solve((RealLinearOperator)AtA, Atb);
		}

		public RealVector calcResidual(RealVector x) {
			return A.operate(x).subtract(b);
		}

		public Errors calcErrors(RealVector x) {
			return new Errors(calcResidual(x));
		}
	}

	public static class Errors {

		public final RealVector residual;
		public final double min;
		public final double max;
		public final double avg;
		public final double rms;

		public Errors(RealVector residual) {

			this.residual = residual;

			double sum = 0.0;
			double sumsq = 0.0;
			double min = Double.POSITIVE_INFINITY;
			double max = Double.NEGATIVE_INFINITY;

			int n = residual.getDimension();
			for (int row=0; row<n; row++) {

				double val = Math.abs(residual.getEntry(row));

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
			return String.format("bounds [%8.4f,%8.4f]   avg %8.4f   rms %8.4f", min, max, avg, rms);
		}
	}

	// TODO: implement triples

	public final SimpleConfSpace confSpace;
	public final EnergyMatrix emat;

	private final Set<RCTuple> tuples = new LinkedHashSet<>();

	public LUTE(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
		this.emat = new EnergyMatrix(confSpace);
	}

	public void addUnprunedPairTuples(PruningMatrix pmat) {

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

						// we found it! It's an unpruned tuple!
						// NOTE: make the tuple in pos2, pos1 order so the positions are already sorted
						// (because pos2 < pos1 by definition)
						addTuple(new RCTuple(pos2, rc2, pos1, rc1));
					}
				}
			}
		}
	}

	public void addTuple(RCTuple tuple) {
		tuple.checkSortedPositions();
		tuples.add(tuple);
	}

	// TODO: make List return
	public ArrayList<RCTuple> getTuples() {
		return new ArrayList<>(tuples);
	}

	public Errors fit(ConfEnergyCalculator confEcalc, ConfDB.ConfTable confTable, int minSamplesPerTuple) {

		// sample the training and test sets
		ConfSampler sampler = new ConfSampler(confSpace, tuples);
		log("sampling confs for %d tuples...", tuples.size());
		Stopwatch samplingSw = new Stopwatch().start();
		Map<RCTuple,Set<int[]>> trainingSet = sampleConfsForTuples(sampler, minSamplesPerTuple);
		Map<RCTuple,Set<int[]>> testSet = sampleConfsForTuples(sampler, minSamplesPerTuple);
		log("done in %s", samplingSw.stop().getTime(2));

		// calculate energies for all the samples
		Map<int[],Double> energies = new Conf.Map<>();
		for (Map<RCTuple,Set<int[]>> set : Arrays.asList(trainingSet, testSet)) {
			for (Set<int[]> samples : set.values()) {
				for (int[] conf : samples) {
					energies.put(conf, null);
				}
			}
		}
		Progress progress = new Progress(energies.size());
		log("calculating energies for %d samples...", progress.getTotalWork());
		for (Map.Entry<int[],Double> entry : energies.entrySet()) {
			int[] conf = entry.getKey();
			confEcalc.calcEnergyAsync(new RCTuple(conf), confTable, (energy) -> {
				entry.setValue(energy);
				progress.incrementProgress();
			});
		}
		confEcalc.tasks.waitForFinish();

		// fit the linear system to the training set
		Stopwatch trainingSw = new Stopwatch().start();
		LinearSystem trainingSystem = new LinearSystem(tuples, trainingSet, energies);
		RealVector x = trainingSystem.fit();
		log("done in %s", trainingSw.stop().getTime(2));

		// analyze the training set errors
		log("training errors: %s", trainingSystem.calcErrors(x));

		// analyze the test set errors
		LinearSystem testSystem = new LinearSystem(tuples, testSet, energies);
		Errors errors = testSystem.calcErrors(x);
		log("test errors:     %s", errors);

		// populate the energy matrix with x
		emat.setConstTerm(0.0);
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
				emat.setOneBody(pos.index, rc.index, 0.0);
			}
		}
		for (int i=0; i<trainingSystem.tuples.size(); i++) {
			RCTuple tuple = trainingSystem.tuples.get(i);
			emat.setPairwise(
				tuple.pos.get(0),
				tuple.RCs.get(0),
				tuple.pos.get(1),
				tuple.RCs.get(1),
				x.getEntry(i)
			);
		}

		return errors;
	}

	private Map<RCTuple,Set<int[]>> sampleConfsForTuples(ConfSampler sampler, int minSamplesPerTuple) {

		// track all samples by tuple
		Map<RCTuple,Set<int[]>> samplesByTuple = new HashMap<>();
		for (RCTuple tuple : tuples) {
			samplesByTuple.put(tuple, new Conf.Set());
		}

		while (true) {

			// find the least sampled tuple
			int minNumSamples = 0;
			RCTuple leastSampledTuple = null;
			for (Map.Entry<RCTuple,Set<int[]>> entry : samplesByTuple.entrySet()) {

				RCTuple tuple = entry.getKey();
				int numSamples = entry.getValue().size();

				if (leastSampledTuple == null || numSamples < minNumSamples) {
					minNumSamples = numSamples;
					leastSampledTuple = tuple;
				}
			}

			// are we done yet?
			int numSamplesNeeded = minSamplesPerTuple - minNumSamples;
			if (numSamplesNeeded <= 0) {
				break;
			}

			// can we get another sample?
			BigInteger maxNumConfs = sampler.getNumConfsUpperBound(leastSampledTuple);
			if (BigInteger.valueOf(minNumSamples).compareTo(maxNumConfs) >= 0) {
				// nope, out of confs to sample, this will have to be good enough
				break;
			}

			// sample another for this tuple (and other tuples too)
			int[] conf = sampler.sample(leastSampledTuple, samplesByTuple.get(leastSampledTuple), 10000);
			assert (conf != null);

			// update the tuple->samples map with the new sample
			for (RCTuple pair : Conf.getPairs(conf)) {
				samplesByTuple.get(pair).add(conf);
			}
		}

		return samplesByTuple;
	}
}
