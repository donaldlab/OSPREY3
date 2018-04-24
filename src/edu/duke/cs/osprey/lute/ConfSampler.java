package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.Conf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TuplesIndex;

import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;


public abstract class ConfSampler {

	public static class Samples {

		public final TuplesIndex tuples;

		private final Set<int[]> confs;
		private final Map<RCTuple,Set<int[]>> confsByTuple;
		private final int[] numConfsByTuple;

		public Samples(TuplesIndex tuples) {

			this.tuples = tuples;

			confs = new Conf.Set();

			confsByTuple = new HashMap<>();
			for (RCTuple tuple : tuples) {
				confsByTuple.put(tuple, new Conf.Set());
			}

			numConfsByTuple = new int[tuples.size()];
			Arrays.fill(numConfsByTuple, 0);
		}

		public RCTuple getLeastSampledTuple(Set<RCTuple> except) {

			// NOTE: using a priority queue to attempt to speed up this query will actually be really slow
			// since we need to update many tuples every time we add a conf
			// and PriorityQueue.remove() is actually linear time =(

			// a simple linear search will have to be fast enough for now
			int bestt = -1;
			int bestSize = -1;
			for (int t=0; t<tuples.size(); t++) {

				int size = numConfsByTuple[t];
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

		// don't let callers edit these sets
		// since we need to keep confsByTuple and confs synchronized

		public Set<int[]> getAllConfs() {
			return Collections.unmodifiableSet(confs);
		}

		public Set<int[]> getConfs(RCTuple tuple) {
			return Collections.unmodifiableSet(confsByTuple.get(tuple));
		}

		public void addConf(int[] conf) {
			tuples.forEachIn(conf, false, true, (index) -> {
				RCTuple tuple = tuples.get(index);
				Set<int[]> confsForTuple = confsByTuple.get(tuple);
				confsForTuple.add(conf);
				numConfsByTuple[index] = confsForTuple.size();
			});
			confs.add(conf);
		}
	}

	public final SimpleConfSpace confSpace;
	public final int randomSeed;

	protected final Random rand;

	public ConfSampler(SimpleConfSpace confSpace, int randomSeed) {

		this.confSpace = confSpace;
		this.randomSeed = randomSeed;

		this.rand = new Random(randomSeed);
	}

	public abstract void sampleConfsForTuples(Samples samples, int minSamplesPerTuple);
}
