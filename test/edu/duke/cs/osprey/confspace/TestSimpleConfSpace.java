package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.minimization.ObjectiveFunction.DofBounds;
import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

public class TestSimpleConfSpace extends TestBase {
	
	private static Molecule mol;
	
	@BeforeClass
	public static void beforeClass() {
		mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
	}
	
	public Strand makeStrand() {
		return Strand.builder(mol)
			.setLovellTemplateLibrary() // explicitly choose Lovell rotamers
			.build();
	}
	
	@Test
	public void moleculeCopy() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		
		// all the residues should be there, but not the same instances
		assertThat(pmol.mol.residues.size(), is(strand.mol.residues.size()));
		for (Residue res : strand.mol.residues) {
			Residue pres = pmol.mol.getResByPDBResNumber(res.getPDBResNumber());
			assertThat(res == pres, is(false));
			assertThat(res.coords == pres.coords, is(false));
		}
		
		// all bonds should be marked, and templates assigned
		for (Residue pres : pmol.mol.residues) {
			assertThat(pres.template, is(not(nullValue())));
			assertThat(pres.intraResBondsMarked, is(true));
			assertThat(pres.interResBondsMarked, is(true));
		}
	}
	
	@Test
	public void onePosition() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		// check the pos
		assertThat(confSpace.positions.size(), is(1));
		SimpleConfSpace.Position pos = confSpace.positions.get(0);
		assertThat(pos.index, is(0));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("2"));
	}
	
	@Test
	public void twoPositions() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		strand.flexibility.get(42).setLibraryRotamers();
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		assertThat(confSpace.positions.size(), is(2));
		
		SimpleConfSpace.Position pos = confSpace.positions.get(0);
		assertThat(pos.index, is(0));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("2"));
		
		pos = confSpace.positions.get(1);
		assertThat(pos.index, is(1));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("42"));
	}
	
	@Test
	public void onePositionShell() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		assertThat(confSpace.shellResNumbers.size(), is(72 - 1));
		
		for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
			assertThat(resNum, not(isIn(confSpace.shellResNumbers)));
		}
		for (String resNum : strand.flexibility.getStaticResidueNumbers()) {
			assertThat(resNum, isIn(confSpace.shellResNumbers));
		}
	}
	
	@Test
	public void twoPositionsShell() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		strand.flexibility.get(42).setLibraryRotamers();
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		assertThat(confSpace.shellResNumbers.size(), is(72 - 2));
		
		for (String resNum : strand.flexibility.getFlexibleResidueNumbers()) {
			assertThat(resNum, not(isIn(confSpace.shellResNumbers)));
		}
		for (String resNum : strand.flexibility.getStaticResidueNumbers()) {
			assertThat(resNum, isIn(confSpace.shellResNumbers));
		}
	}
	
	@Test
	public void oneDiscreteWildType() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers();
		assertThat(strand.mol.getResByPDBResNumber("2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(1));
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("ALA")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(nullValue()));
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
	}
	
	@Test
	public void oneDiscreteGlycine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers("GLY");
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(1));
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("GLY")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(nullValue()));
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("GLY"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
	}
	
	@Test
	public void oneDiscreteWildTypeAndGlycine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers(Strand.WildType, "GLY");
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(2));
		
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("ALA")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(nullValue()));
		
		resConf = resConfs.get(1);
		assertThat(resConf.index, is(1));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("GLY")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(nullValue()));
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("GLY"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
	}
	
	@Test
	public void oneDiscreteValine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers("VAL");
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(3)); // valine has 3 rotamers
		
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("VAL")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(0));
		
		resConf = resConfs.get(1);
		assertThat(resConf.index, is(1));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("VAL")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(1));
		
		resConf = resConfs.get(2);
		assertThat(resConf.index, is(2));
		assertThat(resConf.template, is(strand.templateLib.getTemplate("VAL")));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(2));
		
		ResidueTemplate template = strand.templateLib.getTemplate("VAL");
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0)); // valine has one chi angle, but this is discrete flex
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		double expChi1 = template.getRotamericDihedrals(0, 0);
		double obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		expChi1 = template.getRotamericDihedrals(1, 0);
		obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		expChi1 = template.getRotamericDihedrals(2, 0);
		obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
	}
	
	@Test
	public void oneContinuousValine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers("VAL").setContinuous(3);
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		ResidueTemplate template = strand.templateLib.getTemplate("VAL");
		
		// check DOFs and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		Residue res = pmol.mol.getResByPDBResNumber("2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1)); // valine has one chi angle
		FreeDihedral dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		DofBounds bounds = confSpace.makeBounds(conf);
		assertThat(bounds.size(), is(1));
		double chi1 = template.getRotamericDihedrals(0, 0);
		assertThat(bounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(bounds.getMax(0), isRelatively(chi1 + 3));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		bounds = confSpace.makeBounds(conf);
		assertThat(bounds.size(), is(1));
		chi1 = template.getRotamericDihedrals(1, 0);
		assertThat(bounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(bounds.getMax(0), isRelatively(chi1 + 3));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		bounds = confSpace.makeBounds(conf);
		assertThat(bounds.size(), is(1));
		chi1 = template.getRotamericDihedrals(2, 0);
		assertThat(bounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(bounds.getMax(0), isRelatively(chi1 + 3));
	}
	
	@Test
	public void wildTypeRotamersOnly() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).addWildTypeRotamers();
		assertThat(strand.mol.getResByPDBResNumber("2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(1)); // no alternates in 1CC8
		
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template.name, is("ALA"));
		assertThat(resConf.template == strand.templateLib.getTemplate("ALA"), is(false));
		assertThat(resConf.type, is(ResidueConf.Type.WildType));
		assertThat(resConf.rotamerIndex, is(0));
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
	}
	
	@Test
	public void wildTypeAndLibraryRotamers() {
		
		Strand strand = makeStrand();
		strand.flexibility.get(2).setLibraryRotamers("VAL").addWildTypeRotamers();
		assertThat(strand.mol.getResByPDBResNumber("2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = SimpleConfSpace.build(strand);
		
		List<ResidueConf> resConfs = confSpace.positions.get(0).resConfs;
		assertThat(resConfs.size(), is(4)); // no alternates in 1CC8, 3 rotamers for valine
		
		// library rotamers always get added first
		
		ResidueConf resConf = resConfs.get(0);
		assertThat(resConf.index, is(0));
		assertThat(resConf.template.name, is("VAL"));
		assertThat(resConf.template == strand.templateLib.getTemplate("VAL"), is(true));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(0));
		
		resConf = resConfs.get(1);
		assertThat(resConf.index, is(1));
		assertThat(resConf.template.name, is("VAL"));
		assertThat(resConf.template == strand.templateLib.getTemplate("VAL"), is(true));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(1));
		
		resConf = resConfs.get(2);
		assertThat(resConf.index, is(2));
		assertThat(resConf.template.name, is("VAL"));
		assertThat(resConf.template == strand.templateLib.getTemplate("VAL"), is(true));
		assertThat(resConf.type, is(ResidueConf.Type.Library));
		assertThat(resConf.rotamerIndex, is(2));
		
		resConf = resConfs.get(3);
		assertThat(resConf.index, is(3));
		assertThat(resConf.template.name, is("ALA"));
		assertThat(resConf.template == strand.templateLib.getTemplate("ALA"), is(false));
		assertThat(resConf.type, is(ResidueConf.Type.WildType));
		assertThat(resConf.rotamerIndex, is(0));
		
		// check sequence, DOFs, and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
		
		conf = new int[] { 3 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(confSpace.makeBounds(conf).size(), is(0));
	}
}
