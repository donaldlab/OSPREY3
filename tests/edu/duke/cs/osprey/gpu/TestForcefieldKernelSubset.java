package edu.duke.cs.osprey.gpu;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.BigForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.energy.forcefield.ResPairEnergy;
import edu.duke.cs.osprey.energy.forcefield.SingleResEnergy;
import edu.duke.cs.osprey.gpu.TestForceFieldKernel.Residues;
import edu.duke.cs.osprey.gpu.cuda.GpuStreamPool;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Residue;

public class TestForcefieldKernelSubset extends TestBase {
	
	private static class Forcefields {
		
		public MultiTermEnergyFunction efunc;
		public BigForcefieldEnergy bigff;
		public GpuQueuePool openclQueuePool;
		public GpuForcefieldEnergy gpuffopencl;
		public MultiTermEnergyFunction efuncSub;
		public BigForcefieldEnergy.Subset bigffSub;
		public GpuForcefieldEnergy gpuffSubopencl;
		public GpuStreamPool cudaContextPool;
		public GpuForcefieldEnergy gpuffcuda;
		public GpuForcefieldEnergy gpuffSubcuda;
		
		public void cleanup() {
			gpuffopencl.cleanup();
			openclQueuePool.cleanup();
			gpuffcuda.cleanup();
		}
	}
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	private Forcefields makeForcefields(Residue[] residues)
	throws IOException {
		
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		
		Forcefields ff = new Forcefields();
		
		// make the all pairs energy functions
		ff.efunc = new MultiTermEnergyFunction();
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		for (int pos1=0; pos1<residues.length; pos1++) {
			
			ff.efunc.addTerm(new SingleResEnergy(residues[pos1], ffparams));
			interactions.addResidue(residues[pos1]);
			
			for (int pos2=0; pos2<pos1; pos2++) {
				
				ff.efunc.addTerm(new ResPairEnergy(residues[pos1], residues[pos2], ffparams));
				interactions.addResiduePair(residues[pos1], residues[pos2]);
			}
		}
		
		ff.bigff = new BigForcefieldEnergy(ffparams, interactions);
		ff.openclQueuePool = new GpuQueuePool(1, 1);
		ff.gpuffopencl = new GpuForcefieldEnergy(ffparams, interactions, ff.openclQueuePool);
		ff.cudaContextPool = new GpuStreamPool(1);
		ff.gpuffcuda = new GpuForcefieldEnergy(ffparams, interactions, ff.cudaContextPool);
		
		// make the subset energy functions (first residue against the rest)
		ff.efuncSub = new MultiTermEnergyFunction();
		ff.efuncSub.addTerm(new SingleResEnergy(residues[0], ffparams));
		for (int pos1=1; pos1<residues.length; pos1++) {
			ff.efuncSub.addTerm(new ResPairEnergy(residues[0], residues[pos1], ffparams));
		}
		
		ForcefieldInteractions subsetInteractions = interactions.makeSubsetByResidue(residues[0]);
		ff.bigffSub = ff.bigff.new Subset(subsetInteractions);
		ff.gpuffSubopencl = new GpuForcefieldEnergy(ff.gpuffopencl, subsetInteractions);
		ff.gpuffSubcuda = new GpuForcefieldEnergy(ff.gpuffcuda, subsetInteractions);
		
		return ff;
	}
	
	private void checkEnergies(Residue[] residues, double allPairsEnergy, double subsetEnergy)
	throws IOException {
		
		Forcefields ff = makeForcefields(residues);
		
		// check the all pairs energy functions
		assertThat(ff.efunc.getEnergy(), isRelatively(allPairsEnergy));
		assertThat(ff.bigff.getEnergy(), isRelatively(allPairsEnergy));
		assertThat(ff.gpuffopencl.getEnergy(), isRelatively(allPairsEnergy));
		assertThat(ff.gpuffcuda.getEnergy(), isRelatively(allPairsEnergy));
		
		// check the subset energy functions
		assertThat(ff.efuncSub.getEnergy(), isRelatively(subsetEnergy));
		assertThat(ff.bigffSub.getEnergy(), isRelatively(subsetEnergy));
		assertThat(ff.gpuffSubopencl.getEnergy(), isRelatively(subsetEnergy));
		assertThat(ff.gpuffSubcuda.getEnergy(), isRelatively(subsetEnergy));
		
		ff.cleanup();
	}

	@Test
	public void testSingleGly()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15 };
		checkEnergies(residues, -4.572136255843063, -4.572136255843063);
	}
	
	@Test
	public void testGlyPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly06, r.gly15 };
		checkEnergies(residues, -9.17380398335906, -4.601667727515996);
	}
	
	@Test
	public void testGlySerPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17 };
		checkEnergies(residues, -9.48559560659799, -2.6911081922156552);
	}
	
	@Test
	public void testTrpPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.trp18, r.trp25 };
		checkEnergies(residues, -12.625574526252965, -6.218018599252964);
	}
	
	@Test
	public void test4Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25 };
		checkEnergies(residues, -23.31199205572296, -2.756905624257449);
	}
	
	@Test
	public void test6Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24 };
		checkEnergies(residues, -52.316176530733166, -2.7906943839799343);
	}
	
	@Test
	public void test10Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34 };
		checkEnergies(residues, -93.33337795127768, -2.7991581273906516);
	}
	
	@Test
	public void test14Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48 };
		checkEnergies(residues, -112.44246575304817, -2.799964297346741);
	}
	
	@Test
	public void test24Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = {
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		};
		checkEnergies(residues, -163.74206898485193, -4.612537058185951);
	}
	
	@Test
	public void test6ResiduesMutation()
	throws Exception {
		
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24 };
		
		double expectedWtEnergy = -2.7906943839799343;
		double expectedMutantEnergy = -0.4636208352426065;
		
		Forcefields ff = makeForcefields(residues);
		
		// check the subset energy functions
		assertThat(ff.efuncSub.getEnergy(), isRelatively(expectedWtEnergy));
		assertThat(ff.bigffSub.getEnergy(), isRelatively(expectedWtEnergy));
		assertThat(ff.gpuffSubopencl.getEnergy(), isRelatively(expectedWtEnergy));
		assertThat(ff.gpuffSubcuda.getEnergy(), isRelatively(expectedWtEnergy));
		
		// mutate a residue
		ResidueTypeDOF mutator = new ResidueTypeDOF(r.gly15);
		mutator.mutateTo("VAL");

		assertThat(ff.efuncSub.getEnergy(), isRelatively(expectedMutantEnergy));
		assertThat(ff.bigffSub.getEnergy(), isRelatively(expectedMutantEnergy));
		assertThat(ff.gpuffSubopencl.getEnergy(), isRelatively(expectedMutantEnergy));
		assertThat(ff.gpuffSubcuda.getEnergy(), isRelatively(expectedMutantEnergy));
		
		ff.cleanup();
	}
}
