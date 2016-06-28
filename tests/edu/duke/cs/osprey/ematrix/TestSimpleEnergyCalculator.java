package edu.duke.cs.osprey.ematrix;

import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;

public class TestSimpleEnergyCalculator extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	@Test
	public void testRigid() {
		test(false);
	}
	
	@Test
	public void testContinuous() {
		test(true);
	}
	
	private void test(boolean doMinimize) {
		
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "test/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 8;
		emConfig.addWtRots = true;
		emConfig.doMinimize = doMinimize;
		SearchProblem search = makeSearchProblem(emConfig);
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(
			EnvironmentVars.curEFcnGenerator,
			search.confSpace,
			search.shellResidues
		);
		SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.AllOnSingles;
		
		final double Epsilon = 1e-6;
		
		double exp;
		double obs;
		
    	for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<search.emat.getNumConfAtPos(pos1); rc1++) {

				exp = search.emat.getOneBody(pos1, rc1);
				obs = ecalc.calcSingle(pos1, rc1, dist).getEnergy();
				assertThat(obs, isRelatively(exp, Epsilon));
				
				// pairwise
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<search.emat.getNumConfAtPos(pos2); rc2++) {
						
						exp = search.emat.getPairwise(pos1, rc1, pos2, rc2);
						obs = ecalc.calcPair(pos1, rc1, pos2, rc2, dist).getEnergy();
						assertThat(obs, isRelatively(exp, Epsilon));
					}
				}
			}	
    	}
	}
}
