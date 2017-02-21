package edu.duke.cs.osprey.dof;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

public class TestFreeDihedral extends TestBase {
	
	// NOTE: these are determined from TestDihedrals
	// they're the best we can hope to do without improving Protractor and RotationMatrix
	private static final double EpsilonDegrees = 1e-8;
	
	private static Molecule mol;
	
	@BeforeClass
	public static void beforeClass() {
		initDefaultEnvironment();
		
		mol = new Strand.Builder(PDBIO.readFile("examples/1CC8/1CC8.ss.pdb")).build().mol;
	}
	
	private Residue makePhe() {
		
		// get a phenylalanine
		Residue phe = mol.getResByPDBResNumber("9");
		assertThat(phe.template.name, is("PHE"));
		assertThat(phe.getNumDihedrals(), is(2));
		
		return new Residue(phe);
	}
	
	public FreeDihedral makePheChi1() {
		return new FreeDihedral(makePhe(), 0);
	}
	
	public FreeDihedral makePheChi2() {
		return new FreeDihedral(makePhe(), 1);
	}
	
	@Test
	public void pheChi1InitialVal() {
		FreeDihedral chi1 = makePheChi1();
		assertThat(chi1.getCurVal(), is(chi1.measureDihedralDegrees()));
	}
	
	@Test
	public void pheChi1CoarseCircle() {
		checkAngle(makePheChi1(), -360);
		checkAngle(makePheChi1(), -270);
		checkAngle(makePheChi1(), -180);
		checkAngle(makePheChi1(), -90);
		checkAngle(makePheChi1(), 0);
		checkAngle(makePheChi1(), 90);
		checkAngle(makePheChi1(), 180);
		checkAngle(makePheChi1(), 270);
		checkAngle(makePheChi1(), 360);
	}

	@Test
	public void pheChi1CoarseCircleChained() {
		FreeDihedral chi1 = makePheChi1();
		checkAngle(chi1, -360);
		checkAngle(chi1, -270);
		checkAngle(chi1, -180);
		checkAngle(chi1, -90);
		checkAngle(chi1, 0);
		checkAngle(chi1, 90);
		checkAngle(chi1, 180);
		checkAngle(chi1, 270);
		checkAngle(chi1, 360);
	}
	
	@Test
	public void pheChi1FineCircleChained() {
		FreeDihedral chi1 = makePheChi1();
		int numSamples = 360*2*1024;
		for (int i=0; i<numSamples; i++) {
			double angleDegrees = (double)i*720/numSamples - 360;
			checkAngle(chi1, angleDegrees);
		}
	}
	
	@Test
	public void pheChi2CoarseCircle() {
		checkAngle(makePheChi2(), -360);
		checkAngle(makePheChi2(), -270);
		checkAngle(makePheChi2(), -180);
		checkAngle(makePheChi2(), -90);
		checkAngle(makePheChi2(), 0);
		checkAngle(makePheChi2(), 90);
		checkAngle(makePheChi2(), 180);
		checkAngle(makePheChi2(), 270);
		checkAngle(makePheChi2(), 360);
	}

	@Test
	public void pheChi2CoarseCircleChained() {
		FreeDihedral chi2 = makePheChi2();
		checkAngle(chi2, -360);
		checkAngle(chi2, -270);
		checkAngle(chi2, -180);
		checkAngle(chi2, -90);
		checkAngle(chi2, 0);
		checkAngle(chi2, 90);
		checkAngle(chi2, 180);
		checkAngle(chi2, 270);
		checkAngle(chi2, 360);
	}
	
	@Test
	public void pheChi2FineCircleChained() {
		FreeDihedral chi2 = makePheChi2();
		int numSamples = 360*2*1024;
		for (int i=0; i<numSamples; i++) {
			double angleDegrees = (double)i*720/numSamples - 360;
			checkAngle(chi2, angleDegrees);
		}
	}
	
	@Test
	public void phiChi1Chi2CoarseCircleChained() {
		FreeDihedral chi1 = makePheChi1();
		FreeDihedral chi2 = makePheChi2();
		checkAngles(chi1, -360, chi2,  360);
		checkAngles(chi1, -270, chi2,  270);
		checkAngles(chi1, -180, chi2,  180);
		checkAngles(chi1,  -90, chi2,   90);
		checkAngles(chi1,    0, chi2,    0);
		checkAngles(chi1,   90, chi2,  -90);
		checkAngles(chi1,  180, chi2, -180);
		checkAngles(chi1,  270, chi2, -270);
		checkAngles(chi1,  360, chi2, -360);
	}
	
	private void checkAngle(FreeDihedral dof, double angleDegrees) {
		dof.apply(angleDegrees);
		assertThat(dof.getCurVal(), is(angleDegrees));
		assertThat(dof.measureDihedralDegrees(), isDegrees(Protractor.normalizeDegrees(angleDegrees), EpsilonDegrees));
	}
	
	private void checkAngles(FreeDihedral dof1, double angleDegrees1, FreeDihedral dof2, double angleDegrees2) {
		dof1.apply(angleDegrees1);
		dof2.apply(angleDegrees2);
		assertThat(dof1.getCurVal(), is(angleDegrees1));
		assertThat(dof1.measureDihedralDegrees(), isDegrees(Protractor.normalizeDegrees(angleDegrees1), EpsilonDegrees));
		assertThat(dof2.getCurVal(), is(angleDegrees2));
		assertThat(dof2.measureDihedralDegrees(), isDegrees(Protractor.normalizeDegrees(angleDegrees2), EpsilonDegrees));
	}
}
