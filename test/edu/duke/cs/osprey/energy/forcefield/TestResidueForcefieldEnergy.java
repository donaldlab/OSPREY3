package edu.duke.cs.osprey.energy.forcefield;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.IOException;
import java.util.EnumMap;
import java.util.Map;

import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.TestForcefieldEnergy.Residues;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class TestResidueForcefieldEnergy extends TestBase {
	
	private static enum EfuncType {
		
		AllPairs {
			@Override
			public ResidueForcefieldEnergy makeff(Molecule mol, Residue[] residues, ForcefieldParams ffparams) {
				ResidueInteractions inters = new ResidueInteractions();
				for (int pos1=0; pos1<residues.length; pos1++) {
					inters.addSingle(residues[pos1].getPDBResNumber());
					for (int pos2=0; pos2<pos1; pos2++) {
						inters.addPair(residues[pos1].getPDBResNumber(), residues[pos2].getPDBResNumber());
					}
				}
				return new ResidueForcefieldEnergy(ffparams, inters, mol);
			}
		},
		SingleAndShell {
			@Override
			public ResidueForcefieldEnergy makeff(Molecule mol, Residue[] residues, ForcefieldParams ffparams) {
				ResidueInteractions inters = new ResidueInteractions();
				inters.addSingle(residues[0].getPDBResNumber());
				for (int pos1=1; pos1<residues.length; pos1++) {
					inters.addPair(residues[0].getPDBResNumber(), residues[pos1].getPDBResNumber());
				}
				return new ResidueForcefieldEnergy(ffparams, inters, mol);
			}
		};
		
		public abstract ResidueForcefieldEnergy makeff(Molecule mol, Residue[] residues, ForcefieldParams ffparams);
	}
	
	private void checkEnergies(Residues r, Residue[] residues, double allPairsEnergy, double singleAndShellEnergy)
	throws IOException {
		
		Map<EfuncType,Double> expectedEnergies = new EnumMap<>(EfuncType.class);
		expectedEnergies.put(EfuncType.AllPairs, allPairsEnergy);
		expectedEnergies.put(EfuncType.SingleAndShell, singleAndShellEnergy);
		
		checkEnergies(r, residues, expectedEnergies);
	}
	
	private void checkEnergies(Residues r, Residue[] residues, Map<EfuncType,Double> expectedEnergies) {
		
		ForcefieldParams ffparams = new ForcefieldParams();
		
		for (EfuncType type : EfuncType.values()) {
			ResidueForcefieldEnergy efunc = type.makeff(r.strand.mol, residues, ffparams);
			assertThat(type.toString(), efunc.getEnergy(), isAbsolutely(expectedEnergies.get(type), 1e-10));
		}
	}
	
	@Test
	public void testSingleGly()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15 };
		checkEnergies(r, residues, -4.572136255843063, -4.572136255843063);
	}
	
	@Test
	public void testGlyPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly06, r.gly15 };
		checkEnergies(r, residues, -9.17380398335906, -4.601667727515996);
	}
	
	@Test
	public void testGlySerPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17 };
		checkEnergies(r, residues, -9.48559560659799, -2.6911081922156552);
	}
	
	@Test
	public void testTrpPair()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.trp18, r.trp25 };
		checkEnergies(r, residues, -12.625574526252965, -6.218018599252964);
	}
	
	@Test
	public void test4Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25 };
		checkEnergies(r, residues, -23.31199205572296, -2.756905624257449);
	}
	
	@Test
	public void test6Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24 };
		checkEnergies(r, residues, -52.316176530733166, -2.7906943839799343);
	}
	
	@Test
	public void test10Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34 };
		checkEnergies(r, residues, -93.33337795127768, -2.7991581273906516);
	}
	
	@Test
	public void test14Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48 };
		checkEnergies(r, residues, -112.44246575304817, -2.799964297346741);
	}
	
	@Test
	public void test24Residues()
	throws Exception {
		Residues r = new Residues();
		Residue[] residues = {
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		};
		checkEnergies(r, residues, -163.74206898485193, -4.612537058185951);
	}
	
	@Test
	public void testBrokenProline()
	throws Exception {
		
		Residues r = new Residues();
		Residue[] residues = { r.gly15 };
		
		// mutate to a proline, which will be broken at this pos
		Residue res = r.gly15;
		res.pucker = new ProlinePucker(r.strand.templateLib, res);
		ResidueTypeDOF.switchToTemplate(r.strand.templateLib, res, "PRO");
		
		assertThat(res.confProblems.size(), is(1));
		checkEnergies(r, residues, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
}
