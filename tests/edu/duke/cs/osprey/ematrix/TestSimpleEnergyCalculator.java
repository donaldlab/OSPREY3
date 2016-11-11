package edu.duke.cs.osprey.ematrix;

import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;

public class TestSimpleEnergyCalculator extends TestBase {
	
	private static enum EgenType {
		
		Cpu {
			@Override
			public EnergyFunctionGenerator make() {
				return EnvironmentVars.curEFcnGenerator;
			}
		},
		Gpu {

			@Override
			public EnergyFunctionGenerator make() {
				return new GpuEnergyFunctionGenerator(EnvironmentVars.curEFcnGenerator.ffParams, new GpuQueuePool(1, 1));
			}
		};
		
		public abstract EnergyFunctionGenerator make();
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	@Test
	public void testRigidCpu() {
		test(false, false, EgenType.Cpu);
	}
	
	@Test
	public void testRigidCpuWithMol() {
		test(false, true, EgenType.Cpu);
	}
	
	@Test
	public void testRigidGpu() {
		test(false, false, EgenType.Gpu);
	}
	
	@Test
	public void testRigidGpuWithMol() {
		test(false, true, EgenType.Gpu);
	}
	
	@Test
	public void testContinuousCpu() {
		test(true, false, EgenType.Cpu);
	}
	
	@Test
	public void testContinuousCpuWithMol() {
		test(true, true, EgenType.Cpu);
	}
	
	@Test
	public void testContinuousGpu() {
		test(true, false, EgenType.Gpu);
	}
	
	@Test
	public void testContinuousGpuWithMol() {
		test(true, true, EgenType.Gpu);
	}
	
	private void test(boolean doMinimize, boolean useMolInstance, EgenType egenType) {
		
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "test/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 3;
		emConfig.addWtRots = true;
		emConfig.doMinimize = doMinimize;
		SearchProblem search = makeSearchProblem(emConfig);
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(
			egenType.make(),
			search.confSpace,
			search.shellResidues,
			SimpleEnergyCalculator.ShellDistribution.AllOnSingles
		);
		
		final double Epsilon = 1e-6;
		
		ParameterizedMoleculeCopy pmol = null;
		if (useMolInstance) {
			pmol = new ParameterizedMoleculeCopy(search.confSpace);
		}
		
		double exp;
		double obs;
		
    	for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<search.emat.getNumConfAtPos(pos1); rc1++) {

				// singles
				exp = search.emat.getOneBody(pos1, rc1);
				obs = ecalc.calcSingle(pos1, rc1, pmol).energy;
				assertThat(obs, isRelatively(exp, Epsilon));
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<search.emat.getNumConfAtPos(pos2); rc2++) {
						
						exp = search.emat.getPairwise(pos1, rc1, pos2, rc2);
						obs = ecalc.calcPair(pos1, rc1, pos2, rc2, pmol).energy;
						assertThat(obs, isRelatively(exp, Epsilon));
					}
				}
			}	
    	}
	}
}
