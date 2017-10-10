package edu.duke.cs.osprey.multistatekstar;

import org.junit.Test;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.MSKStarDoer;

public class TestSublinearSearch {

	@Test
    public void test2RL0NegativeDesign() {
		
		String[] args = new String[] {
				"-c", 
				"test/2RL0-neg-des.mskstar/cfgKStar.txt",
				"doMSKStar",
				"test/2RL0-neg-des.mskstar/cfgSystem.txt",
				"test/2RL0-neg-des.mskstar/multistate-spec-sublinear.txt"
	        };
		
		EnvironmentVars.setDataDir("dataFiles");
		
		MSKStarDoer msksd = new MSKStarDoer(args);
		msksd.calcBestSequences();
	}
	
}
