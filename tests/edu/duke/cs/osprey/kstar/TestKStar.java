package edu.duke.cs.osprey.kstar;

import org.junit.Test;

import edu.duke.cs.osprey.control.KStarCalculator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;

public class TestKStar {

	@Test
	public void test2RL0() {
		
		KSConfigFileParser cfp = new KSConfigFileParser(new String[] {
			"-c",
			"test/2RL0.kstar/cfgKStar.txt", 
			"Dummy command",
			"test/2RL0.kstar/cfgMutSearch.txt",
			"test/2RL0.kstar/cfgSystem.txt"
		});
		cfp.loadData();
		
		// configure parallelism
		ThreadParallelism.setNumThreads(cfp.getParams().getInt("numThreads", ThreadParallelism.getNumThreads()));
		MultiTermEnergyFunction.setNumThreads(ThreadParallelism.getNumThreads());

		// run K*
		KStarCalculator kscalc = new KStarCalculator(cfp);
		kscalc.calcKStarScores();
	}
}
