package edu.duke.cs.osprey.energy.forcefield;

import static org.junit.Assert.*;
import static org.hamcrest.Matchers.*;

import java.io.IOException;
import java.util.EnumMap;
import java.util.Map;

import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.dof.ProlinePucker;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

public class TestForcefieldEnergy extends TestBase {
	
	public static class TestResidues {
		
		public final Strand strand;
		public final Residue gly06, gly15, ser17, trp18, trp25, arg22, ala24, ile26, phe31, arg32, glu34,
			val36, leu39, trp47, leu48, ile53, arg55, val56, leu57, ile59, val62, leu64, val65, met66;
		
		public TestResidues() {
			strand = new Strand.Builder(PDBIO.readFile("examples/DAGK/2KDC.P.forOsprey.pdb"))
				.setErrorOnNonTemplateResidues(true)
				.build();
			Molecule mol = strand.mol;
			gly06 = mol.getResByPDBResNumber("6");
			gly15 = mol.getResByPDBResNumber("15");
			ser17 = mol.getResByPDBResNumber("17");
			trp18 = mol.getResByPDBResNumber("18");
			trp25 = mol.getResByPDBResNumber("25");
			arg22 = mol.getResByPDBResNumber("22");
			ala24 = mol.getResByPDBResNumber("24");
			ile26 = mol.getResByPDBResNumber("26");
			phe31 = mol.getResByPDBResNumber("31");
			arg32 = mol.getResByPDBResNumber("32");
			glu34 = mol.getResByPDBResNumber("34");
			val36 = mol.getResByPDBResNumber("36");
			leu39 = mol.getResByPDBResNumber("39");
			trp47 = mol.getResByPDBResNumber("47");
			leu48 = mol.getResByPDBResNumber("48");
			ile53 = mol.getResByPDBResNumber("53");
			arg55 = mol.getResByPDBResNumber("55");
			val56 = mol.getResByPDBResNumber("56");
			leu57 = mol.getResByPDBResNumber("57");
			ile59 = mol.getResByPDBResNumber("59");
			val62 = mol.getResByPDBResNumber("62");
			leu64 = mol.getResByPDBResNumber("64");
			val65 = mol.getResByPDBResNumber("65");
			met66 = mol.getResByPDBResNumber("66");
		}
	}
	
	private static enum EfuncType {
		
		AllPairs {
			@Override
			public MultiTermEnergyFunction makeff(Residue[] residues, ForcefieldParams ffparams) {
				MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
				for (int pos1=0; pos1<residues.length; pos1++) {
					efunc.addTerm(new SingleResEnergy(residues[pos1], ffparams));
					for (int pos2=0; pos2<pos1; pos2++) {
						efunc.addTerm(new ResPairEnergy(residues[pos1], residues[pos2], ffparams));
					}
				}
				return efunc;
			}
		},
		SingleAndShell {
			@Override
			public MultiTermEnergyFunction makeff(Residue[] residues, ForcefieldParams ffparams) {
				MultiTermEnergyFunction efunc = new MultiTermEnergyFunction();
				efunc.addTerm(new SingleResEnergy(residues[0], ffparams));
				for (int pos1=1; pos1<residues.length; pos1++) {
					efunc.addTerm(new ResPairEnergy(residues[0], residues[pos1], ffparams));
				}
				return efunc;
			}
		};
		
		public abstract MultiTermEnergyFunction makeff(Residue[] residues, ForcefieldParams ffparams);
	}
	
	private void checkEnergies(Residue[] residues, double allPairsEnergy, double singleAndShellEnergy)
	throws IOException {
		
		Map<EfuncType,Double> expectedEnergies = new EnumMap<>(EfuncType.class);
		expectedEnergies.put(EfuncType.AllPairs, allPairsEnergy);
		expectedEnergies.put(EfuncType.SingleAndShell, singleAndShellEnergy);
		
		checkEnergies(residues, expectedEnergies);
	}
	
	private void checkEnergies(Residue[] residues, Map<EfuncType,Double> expectedEnergies) {
		
		ForcefieldParams ffparams = new ForcefieldParams();
		
		for (EfuncType type : EfuncType.values()) {
			MultiTermEnergyFunction efunc = type.makeff(residues, ffparams);
			assertThat(type.toString(), efunc.getEnergy(), isAbsolutely(expectedEnergies.get(type)));
		}
	}
	
	@Test
	public void testSingleGly()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15 };
		checkEnergies(residues, -4.572136255843063, -4.572136255843063);
	}
	
	@Test
	public void testGlyPair()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly06, r.gly15 };
		checkEnergies(residues, -9.17380398335906, -4.601667727515996);
	}
	
	@Test
	public void testGlySerPair()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15, r.ser17 };
		checkEnergies(residues, -9.48559560659799, -2.6911081922156552);
	}
	
	@Test
	public void testTrpPair()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.trp18, r.trp25 };
		checkEnergies(residues, -12.625574526252965, -6.218018599252964);
	}
	
	@Test
	public void test4Residues()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25 };
		checkEnergies(residues, -23.31199205572296, -2.756905624257449);
	}
	
	@Test
	public void test6Residues()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24 };
		checkEnergies(residues, -52.316176530733166, -2.7906943839799343);
	}
	
	@Test
	public void test10Residues()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34 };
		checkEnergies(residues, -93.33337795127768, -2.7991581273906516);
	}
	
	@Test
	public void test14Residues()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36, r.leu39, r.trp47, r.leu48 };
		checkEnergies(residues, -112.44246575304817, -2.799964297346741);
	}
	
	@Test
	public void test24Residues()
	throws Exception {
		TestResidues r = new TestResidues();
		Residue[] residues = {
			r.gly06, r.gly15, r.ser17, r.trp18, r.trp25, r.arg22, r.ala24, r.ile26, r.phe31, r.arg32, r.glu34, r.val36,
			r.leu39, r.trp47, r.leu48, r.ile53, r.arg55, r.val56, r.leu57, r.ile59, r.val62, r.leu64, r.val65, r.met66
		};
		checkEnergies(residues, -163.74206898485193, -4.612537058185951);
	}
	
	@Test
	public void testBrokenProline()
	throws Exception {
		
		TestResidues r = new TestResidues();
		Residue[] residues = { r.gly15 };
		
		// mutate to a proline, which will be broken at this pos
		Residue res = r.gly15;
		res.pucker = new ProlinePucker(r.strand.templateLib, res);
		ResidueTypeDOF.switchToTemplate(r.strand.templateLib, res, "PRO");
		
		assertThat(res.confProblems.size(), is(1));
		checkEnergies(residues, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
	}
}
