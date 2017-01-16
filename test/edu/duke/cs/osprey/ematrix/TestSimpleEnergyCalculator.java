package edu.duke.cs.osprey.ematrix;

import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;

public class TestSimpleEnergyCalculator extends TestBase {
	
	public static enum Type {
		
		Cpu {
			
			@Override
			public SimpleEnergyCalculator make(SearchProblem search) {
				return new SimpleEnergyCalculator.Cpu(ffparams, search.confSpace, search.shellResidues);
			}
			
			@Override
			public void cleanup() {
				// nothing to do
			}
			
		},
		Gpu {
			
			private GpuStreamPool pool;
			
			@Override
			public SimpleEnergyCalculator make(SearchProblem search) {
				pool = new GpuStreamPool(1, 1);
				return new SimpleEnergyCalculator.Cuda(pool, ffparams, search.confSpace, search.shellResidues);
			}
			
			@Override
			public void cleanup() {
				pool.cleanup();
			}
		};
		
		public abstract SimpleEnergyCalculator make(SearchProblem search);
		public abstract void cleanup();
	}
	
	private static ForcefieldParams ffparams;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		ffparams = makeDefaultFFParams();
	}
	
	@Test
	public void testRigidCpu() {
		test(false, false, Type.Cpu);
	}
	
	@Test
	public void testRigidCpuWithMol() {
		test(false, true, Type.Cpu);
	}
	
	@Test
	public void testRigidGpu() {
		test(false, false, Type.Gpu);
	}
	
	@Test
	public void testRigidGpuWithMol() {
		test(false, true, Type.Gpu);
	}
	
	@Test
	public void testContinuousCpu() {
		test(true, false, Type.Cpu);
	}
	
	@Test
	public void testContinuousCpuWithMol() {
		test(true, true, Type.Cpu);
	}
	
	@Test
	public void testContinuousGpu() {
		test(true, false, Type.Gpu);
	}
	
	@Test
	public void testContinuousGpuWithMol() {
		test(true, true, Type.Gpu);
	}
	
	private void test(boolean doMinimize, boolean useMolInstance, Type type) {
		
		EnergyMatrixConfig emConfig = new EnergyMatrixConfig();
		emConfig.pdbPath = "examples/DAGK/2KDC.P.forOsprey.pdb";
		emConfig.numFlexible = 3;
		emConfig.addWtRots = true;
		emConfig.doMinimize = doMinimize;
		SearchProblem search = makeSearchProblem(emConfig);
		
		SimpleEnergyCalculator ecalc = type.make(search);
		
		ParameterizedMoleculeCopy pmol;
		if (useMolInstance) {
			pmol = new ParameterizedMoleculeCopy(search.confSpace);
		} else {
			pmol = ParameterizedMoleculeCopy.makeNoCopy(search.confSpace);
		}
		
		final double Epsilon;
		if (doMinimize) {
			Epsilon = 1e-10;
		} else {
			Epsilon = 1e-12;
		}
		
		double exp;
		double efunc;
		double calc;
		
    	for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<search.emat.getNumConfAtPos(pos1); rc1++) {

				// singles
				exp = search.emat.getOneBody(pos1, rc1);
				calc = ecalc.calcSingle(pos1, rc1, pmol).energy;
				assertThat(calc, isAbsolutely(exp, Epsilon));
				
				if (!doMinimize) {
					efunc = ecalc.makeSingleEfunc(pos1, pmol.getCopiedMolecule()).getEnergy();
					assertThat(efunc, isAbsolutely(exp, Epsilon));
				}
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<search.emat.getNumConfAtPos(pos2); rc2++) {
						
						exp = search.emat.getPairwise(pos1, rc1, pos2, rc2);
						calc = ecalc.calcPair(pos1, rc1, pos2, rc2, pmol).energy;
						assertThat(calc, isAbsolutely(exp, Epsilon));
						
						if (!doMinimize) {
							efunc = ecalc.makePairEfunc(pos1, pos2, pmol.getCopiedMolecule()).getEnergy();
							assertThat(efunc, isAbsolutely(exp, Epsilon));
						}
					}
				}
			}	
    	}
    	
    	type.cleanup();
	}
}
