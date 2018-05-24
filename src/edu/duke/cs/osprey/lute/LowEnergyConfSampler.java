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
