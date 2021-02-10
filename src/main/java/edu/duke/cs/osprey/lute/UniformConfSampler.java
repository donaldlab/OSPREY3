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

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TuplesIndex;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

import java.math.BigInteger;
import java.util.*;


/**
 * samples conformations for tuples uniformly from the available possibilities
 * but can get stuck without sampling any confs at all if conformation space is highly pruned
 *
 * don't use this anymore, use the RandomizedDFSConfSampler instead, which solves this problem
 */
@Deprecated
public class UniformConfSampler extends ConfSampler {

	public UniformConfSampler(SimpleConfSpace confSpace, PruningMatrix pmat, int randomSeed) {
		super(confSpace, pmat, randomSeed);
	}

	/**
	 * working with a tuple set (which may not completely cover the conf space)
	 * makes getting an exact count hard, but we can at least get an upper
	 * bound somewhat efficiently
	 */
	public BigInteger getNumConfsUpperBound(RCTuple tuple, TuplesIndex tuples) {

		// small conf space? the best we can do is one conf (the tuple is an entire conf)
		if (tuple.size() == confSpace.positions.size()) {
			return BigInteger.ONE;
		}

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

				if (tuples.isAssignmentCoveredByPairs(conf, pos.index, rc.index)) {
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
	public void sampleConfsForTuples(Samples samples, int minSamplesPerTuple) {

		// keep track of tuples we can't sample anymore
		Set<RCTuple> unsampleableTuples = new HashSet<>();

		while (true) {

			RCTuple tuple = samples.getLeastSampledTuple(unsampleableTuples);
			if (tuple == null) {
				throw new Error("can't find another tuple to sample. This is probably a bug?");
			}
			int numSamples = samples.getConfs(tuple).size();

			// are we done yet?
			if (numSamples >= minSamplesPerTuple) {
				break;
			}

			// can we even get another sample?
			BigInteger maxNumConfs = getNumConfsUpperBound(tuple, samples.tuples);
			if (BigInteger.valueOf(numSamples).compareTo(maxNumConfs) >= 0) {
				// nope, out of confs to sample, that will have to be good enough for this tuple
				unsampleableTuples.add(tuple);
				continue;
			}

			// sample another conf for this tuple (and other tuples too)
			int[] conf = sample(
				tuple,
				samples.getConfs(tuple),
				100,
				samples.tuples
			);
			if (conf == null) {
				// too hard to sample another conf, what we have so far will have to be good enough for this tuple
				// TODO: could maybe do DFS here to sample a conf? not sure how much it will improve fit quality tho
				// since hard sampling means we're probably close to exhausting the conf space
				unsampleableTuples.add(tuple);
				continue;
			}

			samples.addConf(conf);
		}
	}

	public int[] sample(RCTuple tuple, Set<int[]> except, int numAttempts, TuplesIndex tuples) {

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

	public int[] sample(RCTuple tuple, TuplesIndex tuples) {

		// start with a conf with just the tuple assignment
		int[] conf = Conf.make(confSpace, tuple);

		// collect the possible RCs for all unassigned positions
		Set<RCTuple> assignments = new HashSet<>();
		for (SimpleConfSpace.Position pos : confSpace.positions) {

			// skip assigned positions
			if (conf[pos.index] != Conf.Unassigned) {
				continue;
			}

			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
				assignments.add(new RCTuple(pos.index, rc.index));
			}
		}

		while (true) {

			// remove assignments that have been pruned
			assignments.removeIf(assignment ->
				!tuples.isAssignmentCoveredByPairs(
					conf,
					assignment.pos.get(0),
					assignment.RCs.get(0)
				)
				|| isPruned(conf, assignment.pos.get(0), assignment.RCs.get(0))
			);

			// did we run out of possibilities?
			if (assignments.isEmpty()) {
				if (Conf.isCompletelyAssigned(conf)) {
					return conf;
				} else {
					return null;
				}
			}

			// pick a random possibility and assign it to the conf
			RCTuple assignment = removeRandom(assignments);
			assert (conf[assignment.pos.get(0)] == Conf.Unassigned);
			conf[assignment.pos.get(0)] = assignment.RCs.get(0);

			// if the conf completely assigned, we're done
			if (Conf.isCompletelyAssigned(conf)) {
				return conf;
			}

			// otherwise, remove all other assignments for this position and keep going
			assignments.removeIf(a ->
				a.pos.get(0) == assignment.pos.get(0)
			);
		}
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

	private boolean isPruned(int[] conf, int nextPos, int nextRC) {

		// singles and pairs pruning is already accounted for,
		// since we will always use all the pair tuples

		// we just need to check pruning by triple tuples
		for (int pos1=1; pos1<pmat.getNumPos(); pos1++) {

			if (nextPos == pos1 || conf[pos1] == Conf.Unassigned) {
				continue;
			}

			for (int pos2=0; pos2<pos1; pos2++) {

				if (nextPos == pos2 || conf[pos2] == Conf.Unassigned) {
					continue;
				}

				if (pmat.isTriplePruned(pos1, conf[pos1], pos2, conf[pos2], nextPos, nextRC)) {
					return true;
				}
			}
		}
		return false;
	}
}
