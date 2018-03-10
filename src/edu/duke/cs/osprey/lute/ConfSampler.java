package edu.duke.cs.osprey.lute;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.Map;
import java.util.Random;
import java.util.Set;

public abstract class ConfSampler {

	public final SimpleConfSpace confSpace;
	public final int randomSeed;

	protected final Random rand;

	public ConfSampler(SimpleConfSpace confSpace, int randomSeed) {

		this.confSpace = confSpace;
		this.randomSeed = randomSeed;

		this.rand = new Random(randomSeed);
	}

	public abstract Map<RCTuple,Set<int[]>> sampleConfsForTuples(Set<RCTuple> tuples, int minSamplesPerTuple);
}
