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

import java.util.*;
import java.util.function.Function;


public abstract class ConfSampler {

	public static class Samples {

		public final TuplesIndex tuples;

		private final Set<int[]> confs = new Conf.Set();
		private final List<Set<int[]>> confsByTuple = new ArrayList<>();
		private final List<Integer> numConfsByTuple = new ArrayList<>();

		public Samples(TuplesIndex tuples) {
			this.tuples = tuples;
		}

		public RCTuple getLeastSampledTuple(Set<RCTuple> except) {

			// NOTE: using a priority queue to attempt to speed up this query will actually be really slow
			// since we need to update many tuples every time we add a conf
			// and PriorityQueue.remove() is actually linear time =(

			// a simple linear search will have to be fast enough for now
			int bestt = -1;
			int bestSize = -1;
			for (int t=0; t<tuples.size(); t++) {

				int size = numConfsByTuple.get(t);
				if (bestt == -1 || size < bestSize) {

					// skip excepted tuples
					if (except.contains(tuples.get(t))) {
						continue;
					}

					bestt = t;
					bestSize = size;
				}
			}

			if (bestt == -1) {
				return null;
			}
			return tuples.get(bestt);
		}

		public int size() {
			return confs.size();
		}

		public int countTuplesWithLessThan(int numConfs) {
			int count = 0;
			for (int tupleConfs : numConfsByTuple) {
				if (tupleConfs < numConfs) {
					count++;
				}
			}
			return count;
		}

		public int countTuplesWithAtLeast(int numConfs) {
			int count = 0;
			for (int tupleConfs : numConfsByTuple) {
				if (tupleConfs >= numConfs) {
					count++;
				}
			}
			return count;
		}

		// don't let callers edit these sets
		// since we need to keep confsByTuple and confs synchronized

		public Set<int[]> getAllConfs() {
			return Collections.unmodifiableSet(confs);
		}

		public Set<int[]> getConfs(RCTuple tuple) {
			return Collections.unmodifiableSet(confsByTuple.get(tuples.getIndex(tuple)));
		}

		public void addConf(int[] conf) {
			boolean throwIfMissingSingle = conf.length == 1;
			boolean throwIfMissingPair = conf.length > 1;
			tuples.forEachIn(conf, throwIfMissingSingle, throwIfMissingPair, (index) -> {
				Set<int[]> confsForTuple = confsByTuple.get(index);
				confsForTuple.add(conf);
				numConfsByTuple.set(index, confsForTuple.size());
			});
			confs.add(conf);
		}

		public void addTuple(RCTuple tuple) {

			Integer tupleIndex = tuples.getIndex(tuple);
			if (tupleIndex == null) {
				throw new IllegalArgumentException("tuple must be indexed");
			} else if (tupleIndex != tuples.size() - 1) {
				throw new IllegalArgumentException("must be most recently added tuple to the index");
			}
			if (tupleIndex < confsByTuple.size()) {
				throw new IllegalArgumentException("can't add the same tuple more than once: " + tuple);
			}

			Function<int[],Boolean> matches = (conf) -> {
				for (int i=0; i<tuple.size(); i++) {
					int pos = tuple.pos.get(i);
					int rc = tuple.RCs.get(i);
					if (conf[pos] != rc) {
						return false;
					}
				}
				return true;
			};

			Set<int[]> tuplesConfs = new Conf.Set();
			confsByTuple.add(tuplesConfs);

			for (int[] conf : confs) {
				if (matches.apply(conf)) {
					tuplesConfs.add(conf);
				}
			}

			numConfsByTuple.add(tuplesConfs.size());
		}
	}

	public static class NoMoreSamplesException extends RuntimeException {
		private static final long serialVersionUID = 3248938824507606370L;
		NoMoreSamplesException() {
			super("Can't fit LUTE model. No more conformations to sample. Try less pruning?");
		}
	}

	public final SimpleConfSpace confSpace;
	public final PruningMatrix pmat; // just used to avoid sampling pruned sequences when using a sparse tuple basis (e.g. triples)
	public final int randomSeed;

	protected final Random rand;

	public ConfSampler(SimpleConfSpace confSpace, PruningMatrix pmat, int randomSeed) {

		this.confSpace = confSpace;
		this.pmat = pmat;
		this.randomSeed = randomSeed;

		this.rand = new Random(randomSeed);
	}

	public abstract void sampleConfsForTuples(Samples samples, int minSamplesPerTuple);
}
