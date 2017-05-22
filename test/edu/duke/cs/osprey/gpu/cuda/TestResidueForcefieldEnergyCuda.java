package edu.duke.cs.osprey.gpu.cuda;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.Map;

import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ResPairCache;
import edu.duke.cs.osprey.energy.forcefield.TestForcefieldEnergy.TestResidues;
import edu.duke.cs.osprey.gpu.cuda.kernels.ResidueForcefieldEnergyCuda;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.structure.Residues;

public class TestResidueForcefieldEnergyCuda extends TestBase {
	
	// TODO: merge with TestForcefieldKernel?
	
	private static enum EfuncType {
		
		AllPairs {
			@Override
			public ResidueInteractions makeInters(Residues residues) {
				ResidueInteractions inters = new ResidueInteractions();
				for (int pos1=0; pos1<residues.size(); pos1++) {
					inters.addSingle(residues.get(pos1).getPDBResNumber());
					for (int pos2=0; pos2<pos1; pos2++) {
						inters.addPair(residues.get(pos1).getPDBResNumber(), residues.get(pos2).getPDBResNumber());
					}
				}
				return inters;
			}
		},
		SingleAndShell {
			@Override
			public ResidueInteractions makeInters(Residues residues) {
				ResidueInteractions inters = new ResidueInteractions();
				inters.addSingle(residues.get(0).getPDBResNumber());
				for (int pos1=1; pos1<residues.size(); pos1++) {
					inters.addPair(residues.get(0).getPDBResNumber(), residues.get(pos1).getPDBResNumber());
				}
				return inters;
			}
		};
		
		public abstract ResidueInteractions makeInters(Residues residues);
	}
	
	private void checkEnergies(TestResidues r, Residues residues, double allPairsEnergy, double singleAndShellEnergy) {
		
		Map<EfuncType,Double> expectedEnergies = new EnumMap<>(EfuncType.class);
		expectedEnergies.put(EfuncType.AllPairs, allPairsEnergy);
		expectedEnergies.put(EfuncType.SingleAndShell, singleAndShellEnergy);
		
		checkEnergies(r, residues, expectedEnergies);
	}
	
	private void checkEnergies(TestResidues r, Residues residues, Map<EfuncType,Double> expectedEnergies) {
		
		ForcefieldParams ffparams = new ForcefieldParams();
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.setResidues(residues)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		
		GpuStreamPool streams = new GpuStreamPool(1, 1);
		
		try {
			for (EfuncType type : EfuncType.values()) {
				
				ResidueInteractions inters = type.makeInters(residues);
				ResidueForcefieldEnergyCuda efunc = new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, residues);
				double energy = efunc.getEnergy();
				efunc.cleanup();
				
				assertThat(type.toString(), energy, isAbsolutely(expectedEnergies.get(type), 1e-10));
			}
			
		} catch (IOException ex) {
			throw new Error(ex);
		} finally {
			try {
				streams.cleanup();
			} catch (Throwable t) {
				t.printStackTrace(System.err);
			}
		}
	}
	
	@Test
	public void testSingleGly() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15);
		checkEnergies(r, residues, -4.572136255843063, -4.572136255843063);
	}
	
	@Test
	public void testGlyPair() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly06, r.gly15);
		checkEnergies(r, residues, -9.17380398335906, -4.601667727515996);
	}
	
	@Test
	public void testGlySerPair() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17);
		checkEnergies(r, residues, -9.48559560659799, -2.6911081922156552);
	}
	
	@Test
	public void testTrpPair() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.trp18, r.trp25);
		checkEnergies(r, residues, -12.625574526252965, -6.218018599252964);
	}
	
	@Test
	public void test4Residues() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25);
		checkEnergies(r, residues, -23.31199205572296, -2.756905624257449);
	}
	
	@Test
	public void test6Residues() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24);
		checkEnergies(r, residues, -52.316176530733166, -2.7906943839799343);
	}
	
	@Test
	public void test10Residues() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34);
		checkEnergies(r, residues, -93.33337795127768, -2.7991581273906516);
	}
	
	@Test
	public void test14Residues() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48);
		checkEnergies(r, residues, -112.44246575304817, -2.799964297346741);
	}
	
	@Test
	public void test24Residues() {
		TestResidues r = new TestResidues();
		Residues residues = new Residues(
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		);
		checkEnergies(r, residues, -163.74206898485193, -4.612537058185951);
	}
	
	@Test
	public void testBrokenProline() {
		
		TestResidues r = new TestResidues();
		Residues residues = new Residues(r.gly15 );
		
		// mutate to a proline, which will be broken at this pos
		Residue res = r.gly15;
		res.pucker = new ProlinePucker(r.strand.templateLib, res);
		ResidueTypeDOF.switchToTemplate(r.strand.templateLib, res, "PRO");
		
		assertThat(res.confProblems.size(), is(1));
		checkEnergies(r, residues, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
	
	private double calcEnergy(Residues residues, ResidueInteractions inters) {
		
		ForcefieldParams ffparams = new ForcefieldParams();
		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.setResidues(residues)
			.build();
		ResPairCache resPairCache = new ResPairCache(ffparams, connectivity);
		
		
		GpuStreamPool streams = new GpuStreamPool(1, 1);
		ResidueForcefieldEnergyCuda efunc = null;
		try {
			
			efunc = new ResidueForcefieldEnergyCuda(streams, resPairCache, inters, residues);
			return efunc.getEnergy();
			
		} catch (IOException ex) {
			throw new Error(ex);
		} finally {
			try {
				if (efunc != null) {
					efunc.cleanup();
				}
				streams.cleanup();
			} catch (Throwable t) {
				t.printStackTrace(System.err);
			}
		}
	}
	
	@Test
	public void oneIntraWeight() {
		
		TestResidues r = new TestResidues();
		ResidueInteractions inters;
		
		// check base value
		final double baseEnergy = -4.5721362558430645;
		inters = new ResidueInteractions();
		inters.addSingle(r.gly15.getPDBResNumber(), 1, 0);
		assertThat(calcEnergy(new Residues(r.gly15), inters), isAbsolutely(baseEnergy));
		
		// test weight
		for (double weight : Arrays.asList(-0.5, -2.0, -1.0, 0.0, 0.5, 2.0)) {
			inters = new ResidueInteractions();
			inters.addSingle(r.gly15.getPDBResNumber(), weight, 0);
			assertThat("weight: " + weight, calcEnergy(new Residues(r.gly15), inters), isAbsolutely(baseEnergy*weight));
		}
	}
	
	@Test
	public void oneIntraOffset() {
		
		TestResidues r = new TestResidues();
		ResidueInteractions inters;
		
		// check base value
		final double baseEnergy = -4.5721362558430645;
		inters = new ResidueInteractions();
		inters.addSingle(r.gly15.getPDBResNumber(), 1, 0);
		assertThat(calcEnergy(new Residues(r.gly15), inters), isAbsolutely(baseEnergy));
		
		// test weight
		for (double offset : Arrays.asList(-0.5, -2.0, -1.0, 0.0, 0.5, 2.0)) {
			inters = new ResidueInteractions();
			inters.addSingle(r.gly15.getPDBResNumber(), 1, offset);
			assertThat("offset: " + offset, calcEnergy(new Residues(r.gly15), inters), isAbsolutely(baseEnergy + offset));
		}
	}
}
