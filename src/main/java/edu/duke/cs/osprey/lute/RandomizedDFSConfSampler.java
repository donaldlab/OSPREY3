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


import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.util.*;


/**
 * Samples confs for a tuple by doing DFS search.
 * But tree traversal order is randomized so samples confs aren't biased
 * by any particular tree traversal order.
 */
public class RandomizedDFSConfSampler extends ConfSampler {

	public RandomizedDFSConfSampler(SimpleConfSpace confSpace, PruningMatrix pmat, int randomSeed) {
		super(confSpace, pmat, randomSeed);
	}

	@Override
	public void sampleConfsForTuples(Samples samples, int minSamplesPerTuple) {

		// keep track of tuples we can't sample anymore
		Set<RCTuple> unsampleableTuples = new HashSet<>();
		boolean sampledSomething = false;

		while (true) {

			// get the next tuple to sample
			RCTuple tuple = samples.getLeastSampledTuple(unsampleableTuples);
			if (tuple == null) {
				// we can't keep going, there aren't any tuples left to sample
				if (sampledSomething) {
					// at least we made some progress, so that's good
					return;
				} else {
					// no progress, this is bad =(
					throw new ConfSampler.NoMoreSamplesException();
				}
			}
			Set<int[]> confs = samples.getConfs(tuple);

			// are we done yet?
			if (confs.size() >= minSamplesPerTuple) {
				break;
			}

			// get the next sample
			int[] conf = randomizedDFS(tuple, confs);
			if (conf == null) {
				unsampleableTuples.add(tuple);
			} else {
				samples.addConf(conf);
				sampledSomething = true;
			}
		}
	}

	private int[] randomizedDFS(RCTuple tuple, Set<int[]> except) {

		// collect all the positions: assigned first, then unassigned
		List<SimpleConfSpace.Position> positions = new ArrayList<>();
		for (int pos : tuple.pos) {
			positions.add(confSpace.positions.get(pos));
		}
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			if (!positions.contains(pos)) {
				positions.add(pos);
			}
		}

		// pick a random permutation for the unassigned positions
		Collections.shuffle(
			positions.subList(tuple.size(), positions.size()),
			rand
		);

		// do DFS
		Stack<ConfIndex> stack = new Stack<>();

		// start with the DFS node representing the tuple
		ConfIndex root = new ConfIndex(positions.size());
		root.numDefined = tuple.size();
		for (int i=0; i<tuple.size(); i++) {
			root.definedPos[i] = tuple.pos.get(i);
			root.definedRCs[i] = tuple.RCs.get(i);
		}
		root.sortDefined();
		root.updateUndefined();
		stack.push(root);

		while (!stack.isEmpty()) {

			ConfIndex node = stack.pop();

			// hit a leaf node? we're done here
			if (node.numDefined == positions.size()) {
				int[] conf = Conf.make(node);
				if (except.contains(conf)) {
					continue;
				} else {
					return conf;
				}
			}

			// pick the next position to expand (randomly)
			SimpleConfSpace.Position pos = positions.get(node.numDefined);

			// expand the RCs (in random order)
			List<SimpleConfSpace.ResidueConf> rcs = new ArrayList<>(pos.resConfs);
			Collections.shuffle(rcs, rand);
			for (SimpleConfSpace.ResidueConf rc : rcs) {

				// if this child was pruned by the pruning matrix, then skip it
				if (isPruned(node, pos.index, rc.index)) {
					continue;
				}

				// otherwise, expand it
				stack.push(node.assign(pos.index, rc.index));
			}
		}

		return null;
	}

	private boolean isPruned(ConfIndex confIndex, int nextPos, int nextRc) {

		// check pairs
		for (int i=0; i<confIndex.numDefined; i++) {
			int pos = confIndex.definedPos[i];
			int rc = confIndex.definedRCs[i];
			assert (pos != nextPos || rc != nextRc);
			if (pmat.getPairwise(pos, rc, nextPos, nextRc)) {
				return true;
			}
		}

		// check triples
		if (pmat.hasHigherOrderTuples()) {

			RCTuple tuple = new RCTuple(0, 0, 0, 0, 0, 0);

			for (int i1=0; i1<confIndex.numDefined; i1++) {
				int pos1 = confIndex.definedPos[i1];
				int rc1 = confIndex.definedRCs[i1];
				assert (pos1 != nextPos || rc1 != nextRc);

				for (int i2=0; i2<i1; i2++) {
					int pos2 = confIndex.definedPos[i2];
					int rc2 = confIndex.definedRCs[i2];
					assert (pos2 != nextPos || rc2 != nextRc);

					tuple.set(pos1, rc1, pos2, rc2, nextPos, nextRc);
					tuple.sortPositions();

					if (pmat.getTuple(tuple)) {
						return true;
					}
				}
			}
		}

		return false;
	}
}
