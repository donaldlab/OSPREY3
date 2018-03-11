package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;

import java.util.*;
import java.util.function.Function;

import static edu.duke.cs.osprey.tools.Log.log;


public class LowEnergyConfSampler extends ConfSampler {

	public final PruningMatrix pmat;
	public final Function<RCs,ConfSearch> astarFactory;

	public LowEnergyConfSampler(SimpleConfSpace confSpace, int randomSeed, PruningMatrix pmat, Function<RCs,ConfSearch> astarFactory) {
		super(confSpace, randomSeed);

		this.pmat = pmat;
		this.astarFactory = astarFactory;
	}

	@Override
	public Map<RCTuple,Set<int[]>> sampleConfsForTuples(Set<RCTuple> tuples, int minSamplesPerTuple) {

		// track all samples by tuple
		Map<RCTuple,Set<int[]>> samplesByTuple = new HashMap<>();
		for (RCTuple tuple : tuples) {
			samplesByTuple.put(tuple, new Conf.Set());
		}

		// keep track of tuples we can't sample anymore
		Set<RCTuple> unsampleableTuples = new HashSet<>();

		while (true) {

			// find the least sampled tuple
			// TODO: use a priority queue here?
			Set<int[]> alreadySampled = null;
			RCTuple leastSampledTuple = null;
			for (Map.Entry<RCTuple,Set<int[]>> entry : samplesByTuple.entrySet()) {

				RCTuple tuple = entry.getKey();
				Set<int[]> tupleSamples = entry.getValue();

				// skip tuples we can't sample
				if (unsampleableTuples.contains(tuple)) {
					continue;
				}

				if (leastSampledTuple == null || tupleSamples.size() < alreadySampled.size()) {
					alreadySampled = tupleSamples;
					leastSampledTuple = tuple;
				}
			}
			assert (leastSampledTuple != null);
			assert (alreadySampled != null);

			// are we done sampling yet?
			int numSamplesNeeded = minSamplesPerTuple - alreadySampled.size();
			if (numSamplesNeeded <= 0) {
				break;
			}

			// sample more confs for this tuple (and other tuples too)
			int[] leastSampledTupleConf = Conf.make(confSpace, leastSampledTuple);
			RCs rcs = new RCs(pmat);
			rcs = new RCs(rcs, (pos, rc) -> {

				// allow all rcs at unassigned positions
				if (leastSampledTupleConf[pos] == Conf.Unassigned) {
					return true;
				}

				// allow only the assigned RC at assigned positions
				return rc == leastSampledTupleConf[pos];
			});

			if (MathTools.isZero(rcs.getNumConformations())) {
				// no confs in the sub-space induced by this tuple
				// that seems really bad
				throw new Error("tuple " + leastSampledTuple + " has no compatible conformations");
			}

			// get a pool of low-energy confs
			final int maxPoolSize = numSamplesNeeded*10;
			ConfSearch astar = astarFactory.apply(rcs);
			List<int[]> confPool = new ArrayList<>();
			for (int i=0; i<maxPoolSize; i++) {
				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					break;
				}
				confPool.add(conf.getAssignments());
			}

			// sample the desired number of samples
			List<int[]> samples = new ArrayList<>();
			for (int i=0; i<confPool.size(); i++) {
				int[] conf = confPool.get(rand.nextInt(confPool.size()));
				boolean isUnique = alreadySampled.add(conf);
				if (isUnique) {
					samples.add(conf);
					if (samples.size() >= numSamplesNeeded) {
						break;
					}
				}
			}

			if (samples.isEmpty()) {
				// ran out of confs to sample
				// TODO: make the pool bigger?
				// TEMP: pretend we can't sample confs for this tuple anymore
				log("tuple %s ran out of samples, have %d, need %d, pool size %d",
					leastSampledTuple, alreadySampled.size(), minSamplesPerTuple, confPool.size()
				);
				unsampleableTuples.add(leastSampledTuple);
				continue;
			}

			// update the tuple->samples map with the new sample
			for (RCTuple tuple : tuples) {
				for (int[] conf : samples) {
					if (Conf.containsTuple(conf, tuple)) {
						samplesByTuple.get(tuple).add(conf);
					}
				}
			}
		}

		return samplesByTuple;
	}
}
