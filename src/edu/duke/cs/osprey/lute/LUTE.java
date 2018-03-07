package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;


public class LUTE {

	// TODO: implement triples

	public final SimpleConfSpace confSpace;

	private final Set<RCTuple> tuples = new LinkedHashSet<>();

	public LUTE(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
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

	public double fit(ConfEnergyCalculator confEcalc) {

		// track all samples by tuple
		Map<RCTuple,Set<int[]>> samplesByTuple = new HashMap<>();
		for (RCTuple tuple : tuples) {
			samplesByTuple.put(tuple, new Conf.Set());
		}

		ConfSampler sampler = new ConfSampler(confSpace, tuples);

		final int minSamplesPerTuple = 10;

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
			// TEMP
			log("least sampled tuple: %s -> %d", leastSampledTuple, minNumSamples);

			if (minNumSamples >= minSamplesPerTuple) {
				break;
			}

			// sample more confs for this tuple
			Set<int[]> confs = sampler.sample(leastSampledTuple, 10, 10000);
			log("tuple: %s   sampled confs: %d/%d  %s", leastSampledTuple, confs.size(), sampler.getNumConfsUpperBound(leastSampledTuple), confs);

			// update the tuple->samples map with the new samples
			for (int[] conf : confs) {
				for (RCTuple pair : Conf.getPairs(conf)) {
					samplesByTuple.get(pair).add(conf);
				}
			}
		}

		// collect all the sampled confs and calculate their energies
		Map<int[],Double> confs = new Conf.Map<>();
		for (Set<int[]> samples : samplesByTuple

		/*
		final double pruningInterval = 0.0;
		final LUTESettings luteSettings = LUTESettings.defaultLUTE();

		TupleExpander expander = new NewConfETupleExpander(confSpace, pruningInterval, luteSettings, confEcalc, pmat);

		setupSamples(tuplesToFit);//set up the training set (in the process, prune tuples that don't provide reasonable energies)

		fitLeastSquares();

		trainingSamples.updateFitVals(fof);
		CVSamples.updateFitVals(fof);

		System.out.println("TRAINING SAMPLES: ");
		trainingSamples.printResids();

		System.out.println("CV SAMPLES: ");
		CVSamples.printResids();

		return CVSamples.totalResid;
		*/

		// TEMP
		return Double.NaN;
	}

	public EnergyMatrix makeEnergyMatrix() {
		// TODO
		return null;
	}
}
