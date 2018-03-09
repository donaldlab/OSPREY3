package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

import java.math.BigInteger;
import java.util.*;


public class ConfSampler {

	public final SimpleConfSpace confSpace;
	public final Set<RCTuple> tuples;

	private final Random rand = new Random(12345); // deterministic random, rather than stochastic

	public ConfSampler(SimpleConfSpace confSpace, Set<RCTuple> tuples) {
		this.confSpace = confSpace;
		this.tuples = tuples;
	}

	/**
	 * working with a tuple set (which may not completely cover the conf space)
	 * makes getting an exact count hard, but we can at least get an upper
	 * bound somewhat efficiently
	 */
	public BigInteger getNumConfsUpperBound(RCTuple tuple) {

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

				if (isCompatiblePairs(conf, new RCTuple(pos.index, rc.index))) {
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

	public Set<int[]> sample(RCTuple tuple, int numSamples, int numAttempts) {

		// don't know how to prune possible assignments based on a list of confs
		// so I don't think we can do better than sample-and-reject here =(
		Set<int[]> confs = new Conf.Set();
		for (int i=0; i<numAttempts; i++) {
			confs.add(sample(tuple));
			if (confs.size() >= numSamples) {
				break;
			}
		}
		return confs;
	}

	public int[] sample(RCTuple tuple, Set<int[]> except, int numAttempts) {

		// don't know how to prune possible assignments based on a list of confs
		// so I don't think we can do better than sample-and-reject here =(
		for (int i=0; i<numAttempts; i++) {
			int[] conf = sample(tuple);
			if (except == null || !except.contains(conf)) {
				return conf;
			}
		}

		return null;
	}

	public int[] sample(RCTuple tuple) {

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
			possibleAssignments.removeIf(possibleAssignment -> !isCompatiblePairs(conf, possibleAssignment));

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

	private boolean isCompatiblePairs(int[] conf, RCTuple possibleAssignment) {

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
