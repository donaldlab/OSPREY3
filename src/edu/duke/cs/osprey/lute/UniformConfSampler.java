package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

import java.math.BigInteger;
import java.util.*;


public class UniformConfSampler extends ConfSampler {

	public UniformConfSampler(SimpleConfSpace confSpace, int randomSeed) {
		super(confSpace, randomSeed);
	}

	/**
	 * working with a tuple set (which may not completely cover the conf space)
	 * makes getting an exact count hard, but we can at least get an upper
	 * bound somewhat efficiently
	 */
	public BigInteger getNumConfsUpperBound(RCTuple tuple, Set<RCTuple> tuples) {

		// make a conf with the tuple assignment
		int[] conf = Conf.make(confSpace, tuple);

		BigInteger numConfs = BigInteger.ZERO;

		for (SimpleConfSpace.Position pos : confSpace.positions) {

			// skip tuple positions
			if (conf[pos.index] != Conf.Unassigned) {
				continue;
			}

			// count the number of possible assignments that are compatible with this tuple at this pos
			int numAssignments = 0;
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

				if (isCompatiblePairs(conf, new RCTuple(pos.index, rc.index), tuples)) {
					numAssignments++;
				}
			}

			if (MathTools.isZero(numConfs)) {
				numConfs = BigInteger.valueOf(numAssignments);
			} else {
				numConfs = numConfs.multiply(BigInteger.valueOf(numAssignments));
			}
		}

		return numConfs;
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
			int minNumSamples = 0;
			RCTuple leastSampledTuple = null;
			for (Map.Entry<RCTuple,Set<int[]>> entry : samplesByTuple.entrySet()) {

				RCTuple tuple = entry.getKey();
				int numSamples = entry.getValue().size();

				// skip tuples we can't sample
				if (unsampleableTuples.contains(tuple)) {
					continue;
				}

				if (leastSampledTuple == null || numSamples < minNumSamples) {
					minNumSamples = numSamples;
					leastSampledTuple = tuple;
				}
			}

			// are we done yet?
			if (minNumSamples >= minSamplesPerTuple) {
				break;
			}

			// can we even get another sample?
			BigInteger maxNumConfs = getNumConfsUpperBound(leastSampledTuple, tuples);
			if (BigInteger.valueOf(minNumSamples).compareTo(maxNumConfs) >= 0) {
				// nope, out of confs to sample, that will have to be good enough for this tuple
				unsampleableTuples.add(leastSampledTuple);
				continue;
			}

			// sample another conf for this tuple (and other tuples too)
			int[] conf = sample(
				leastSampledTuple,
				samplesByTuple.get(leastSampledTuple),
				minSamplesPerTuple*100,
				tuples
			);
			if (conf == null) {
				// too hard to sample another conf, what we have so far will have to be good enough for this tuple
				// TODO: could maybe do DFS here to sample a conf? not sure how much it will improve fit quality tho
				// since hard sampling means we're probably close to exhausting the conf space
				unsampleableTuples.add(leastSampledTuple);
				continue;
			}

			// update the tuple->samples map with the new sample
			for (RCTuple tuple : tuples) {
				if (Conf.containsTuple(conf, tuple)) {
					samplesByTuple.get(tuple).add(conf);
				}
			}
		}

		return samplesByTuple;
	}

	public Set<int[]> sample(RCTuple tuple, int numSamples, int numAttempts, Set<RCTuple> tuples) {

		// don't know how to prune possible assignments based on a list of confs
		// so I don't think we can do better than sample-and-reject here =(
		Set<int[]> confs = new Conf.Set();
		for (int i=0; i<numAttempts; i++) {
			confs.add(sample(tuple, tuples));
			if (confs.size() >= numSamples) {
				break;
			}
		}
		return confs;
	}

	public int[] sample(RCTuple tuple, Set<int[]> except, int numAttempts, Set<RCTuple> tuples) {

		// don't know how to prune possible assignments based on a list of confs
		// so I don't think we can do better than sample-and-reject here =(
		for (int i=0; i<numAttempts; i++) {
			int[] conf = sample(tuple, tuples);
			if (conf == null) {
				// must have sampled a dead-end
				continue;
			}
			if (except == null || !except.contains(conf)) {
				return conf;
			}
		}

		return null;
	}

	public int[] sample(RCTuple tuple, Set<RCTuple> tuples) {

		// start with a conf with just the tuple assignment
		int[] conf = Conf.make(confSpace, tuple);

		// collect the possible RCs for all unassigned positions
		Set<RCTuple> possibleAssignments = new HashSet<>();
		for (SimpleConfSpace.Position pos : confSpace.positions) {

			// skip assigned positions
			if (conf[pos.index] != Conf.Unassigned) {
				continue;
			}

			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
				possibleAssignments.add(new RCTuple(pos.index, rc.index));
			}
		}

		while (true) {

			// keep only RCs compatible with the conf
			possibleAssignments.removeIf(possibleAssignment -> !isCompatiblePairs(conf, possibleAssignment, tuples));

			// did we run out of possibilities?
			if (possibleAssignments.isEmpty()) {
				if (Conf.isCompletelyAssigned(conf)) {
					return conf;
				} else {
					return null;
				}
			}

			// pick a random possibility and assign it to the conf
			RCTuple assignment = removeRandom(possibleAssignments);
			assert (conf[assignment.pos.get(0)] == Conf.Unassigned);
			conf[assignment.pos.get(0)] = assignment.RCs.get(0);

			// if the conf completely assigned, we're done
			if (Conf.isCompletelyAssigned(conf)) {
				return conf;
			}
		}
	}

	private boolean isCompatiblePairs(int[] conf, RCTuple possibleAssignment, Set<RCTuple> tuples) {

		// all pairs between this possible assignment and the conf assignments must be present in the tuple set
		for (SimpleConfSpace.Position pos : confSpace.positions) {

			// skip unassigned positions
			if (conf[pos.index] == Conf.Unassigned) {
				continue;
			}

			// is this tuple present in the set?
			RCTuple tuple = possibleAssignment.addRC(pos.index, conf[pos.index]).sorted();
			if (!tuples.contains(tuple)) {
				return false;
			}
		}

		return true;
	}

	private <T> T removeRandom(Set<T> set) {

		// random access in an unordered collection in Java is sadly linear time =(

		// since we can't do better than linear, use the iterator to impose an
		// order on the set and remove the item at a random position in the order

		int removeIndex = rand.nextInt(set.size());

		Iterator<T> iter = set.iterator();
		int index = 0;
		while (iter.hasNext()) {
			T item = iter.next();
			if (index == removeIndex) {
				iter.remove();
				return item;
			}
			index++;
		}

		throw new UnpossibleError();
	}
}
