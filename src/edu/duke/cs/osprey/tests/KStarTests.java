package edu.duke.cs.osprey.tests;

import edu.duke.cs.osprey.control.KStarCalculator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.kstar.KSConfigFileParser;
import edu.duke.cs.osprey.parallelism.ThreadParallelism;
import junit.framework.TestCase;
import junit.framework.TestResult;

public class KStarTests extends TestCase {

	KStarCalculator ksc;
	KSConfigFileParser cfp;

	public static void main(String[] args) {
		try {
			
			KStarTests tks = new KStarTests();
			tks.setUp();
			tks.run();
			
		} catch(Exception e) {
			
			System.out.println(e.getMessage());
			e.printStackTrace();
			System.exit(1);
			
		}
	}

	protected void setUp() throws Exception {
		
		super.setUp();
		
		String[] testArgs = new String[]{"-c", "test/2RL0.kstar/cfgKStar.txt", 
				"Dummy command", "test/2RL0.kstar/cfgMutSearch.txt", "test/2RL0.kstar/cfgSystem.txt"};
		cfp = new KSConfigFileParser(testArgs);//args 1, 3+ are configuration files
		cfp.loadData();
		
		ThreadParallelism.setNumThreads(cfp.getParams().getInt("numThreads", ThreadParallelism.getNumThreads()));
		MultiTermEnergyFunction.setNumThreads(ThreadParallelism.getNumThreads());

		ksc = new KStarCalculator(cfp);
		
	}

	public TestResult run() {
		
		ksc.calcKStarScores();
		return null;
		
	}
}
