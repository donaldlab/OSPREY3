/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

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

	public final Function<RCs,ConfSearch> astarFactory;

	public LowEnergyConfSampler(SimpleConfSpace confSpace, PruningMatrix pmat, int randomSeed, Function<RCs,ConfSearch> astarFactory) {
		super(confSpace, pmat, randomSeed);

		this.astarFactory = astarFactory;
	}

	@Override
	public void sampleConfsForTuples(Samples samples, int minSamplesPerTuple) {

		// keep track of tuples we can't sample anymore
		Set<RCTuple> unsampleableTuples = new HashSet<>();

		while (true) {

			RCTuple tuple = samples.getLeastSampledTuple(unsampleableTuples);
			Set<int[]> alreadySampled = samples.getConfs(tuple);

			// are we done sampling yet?
			int numSamplesNeeded = minSamplesPerTuple - alreadySampled.size();
			if (numSamplesNeeded <= 0) {
				break;
			}

			// sample more confs for this tuple (and other tuples too)
			int[] sampleConf = Conf.make(confSpace, tuple);
			RCs rcs = new RCs(pmat);
			rcs = new RCs(rcs, (pos, rc) -> {

				// allow all rcs at unassigned positions
				if (sampleConf[pos] == Conf.Unassigned) {
					return true;
				}

				// allow only the assigned RC at assigned positions
				return rc == sampleConf[pos];
			});

			if (MathTools.isZero(rcs.getNumConformations())) {
				// no confs in the sub-space induced by this tuple
				// that seems really bad
				throw new Error("tuple " + tuple + " has no compatible conformations");
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

			// try to sample the desired number of samples
			int numAttempts = numSamplesNeeded*10;
			for (int i=0; i<numAttempts; i++) {
				int[] conf = confPool.get(rand.nextInt(confPool.size()));
				samples.addConf(conf);
				if (alreadySampled.size() >= numSamplesNeeded) {
					break;
				}
			}

			if (alreadySampled.size() < numSamplesNeeded) {
				// ran out of confs to sample
				// TODO: make the pool bigger?
				// TEMP: pretend we can't sample confs for this tuple anymore
				log("tuple %s ran out of samples, have %d, need %d, pool size %d",
					tuple, alreadySampled.size(), minSamplesPerTuple, confPool.size()
				);
				unsampleableTuples.add(tuple);
			}
		}
	}
}
