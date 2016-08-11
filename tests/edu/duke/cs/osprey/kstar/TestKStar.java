package edu.duke.cs.osprey.kstar;

import org.junit.Test;

import edu.duke.cs.osprey.control.KStarCalculator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import junit.framework.TestResult;

public class TestKStar {

	KStarCalculator ksc;
	KSConfigFileParser cfp;

	@Test
	public void test2RL0() {
		try {
			
			TestKStar tks = new TestKStar();
			tks.setUp();
			tks.run();
			
		} catch(Exception e) {
			
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
			
		}
	}

	protected void setUp() throws Exception {
		
		String[] testArgs = new String[]{"-c", "test/2RL0.kstar/cfgKStar.txt", 
				"Dummy command", "test/2RL0.kstar/cfgMutSearch.txt", "test/2RL0.kstar/cfgSystem.txt"};
		cfp = new KSConfigFileParser(testArgs);//args 1, 3+ are configuration files
		cfp.loadData();
		cfp.getParams().setValue("KSTARPFUNCMETHOD", "parallel2");
		cfp.getParams().setValue("KSTARPFUNCTHREADS", "3");
		cfp.getParams().setValue("NUMTHREADS", "3");
		
		ThreadParallelism.setNumThreads(cfp.getParams().getInt("numThreads", ThreadParallelism.getNumThreads()));
		MultiTermEnergyFunction.setNumThreads(ThreadParallelism.getNumThreads());

		ksc = new KStarCalculator(cfp);
		
	}

	public TestResult run() {
		
		ksc.calcKStarScores();
		return null;
		
	}
}
