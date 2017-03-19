package edu.duke.cs.osprey.gpu;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gpu.TestForceFieldKernel.EnergyFunctionType;
import edu.duke.cs.osprey.gpu.TestForceFieldKernel.Forcefields;
import edu.duke.cs.osprey.gpu.TestForceFieldKernel.Residues;
import edu.duke.cs.osprey.structure.Residue;

public class TestForcefieldParams extends TestBase {
	
	private static final double Epsilon = 1e-12;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
	}
	
	private void checkEnergies(Residue[] residues, double noSolvEnergy, double eef1Energy)
	throws IOException {
		
		ForcefieldParams ffparams = makeDefaultFFParams();
		
		// no solv
		ffparams.doSolvationE = false;
		Forcefields ff = TestForceFieldKernel.makeForcefields(residues, EnergyFunctionType.AllPairs, ffparams);
		assertThat(ff.efunc.getEnergy(), isAbsolutely(noSolvEnergy, Epsilon));
		assertThat(ff.bigff.getEnergy(), isAbsolutely(noSolvEnergy, Epsilon));
		assertThat(ff.gpuffopencl.getEnergy(), isAbsolutely(noSolvEnergy, Epsilon));
		assertThat(ff.gpuffcuda.getEnergy(), isAbsolutely(noSolvEnergy, Epsilon));
		ff.cleanup();
		
		// EEF1 solv
		ffparams.doSolvationE = true;
		ff = TestForceFieldKernel.makeForcefields(residues, EnergyFunctionType.AllPairs, ffparams);
		assertThat(ff.efunc.getEnergy(), isAbsolutely(eef1Energy, Epsilon));
		assertThat(ff.bigff.getEnergy(), isAbsolutely(eef1Energy, Epsilon));
		assertThat(ff.gpuffopencl.getEnergy(), isAbsolutely(eef1Energy, Epsilon));
		assertThat(ff.gpuffcuda.getEnergy(), isAbsolutely(eef1Energy, Epsilon));
		ff.cleanup();
	}
	
	@Test
	public void testSingleGly()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15 };
		checkEnergies(residues, 0.7914065324002029, -4.572136255843063);
	}
	
	@Test
	public void testGlyPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly06, r.gly15 };
		checkEnergies(residues, 1.5576201272528598, -9.17380398335906);
	}
	
	@Test
	public void testGlySerPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17 };
		checkEnergies(residues, 3.026413205967965, -9.48559560659799);
	}
	
	@Test
	public void testTrpPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.trp18, r.trp25 };
		checkEnergies(residues, 3.6482861076947595, -12.625574526252965);
	}
	
	@Test
	public void test4Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25 };
		checkEnergies(residues, 2.970182274508187, -23.31199205572296);
	}
	
	@Test
	public void test6Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24 };
		checkEnergies(residues, -9.856849417475418, -52.316176530733166);
	}
	
	@Test
	public void test10Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34 };
		checkEnergies(residues, -17.302367217342816, -93.33337795127768);
	}
	
	@Test
	public void test14Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48 };
		checkEnergies(residues, -21.693144259934773, -112.44246575304817);
	}
	
	@Test
	public void test24Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = {
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		};
		checkEnergies(residues, -34.032360664731826, -163.74206898485193);
	}
}
