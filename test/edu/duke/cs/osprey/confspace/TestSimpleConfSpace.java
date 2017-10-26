package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.util.List;

import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;
import edu.duke.cs.osprey.dof.FreeDihedral;
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
		return new Strand.Builder(mol)
			// explicitly choose Lovell rotamers
			.setTemplateLibrary(new ResidueTemplateLibrary.Builder()
				.clearRotamers()
				.addLovellRotamers()
				.build())
			.build();
	}
	
	@Test
	public void moleculeCopy() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
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
		
		// and the indices should be set
		for (int i=0; i<pmol.mol.residues.size(); i++) {
			assertThat(pmol.mol.residues.get(i).indexInMolecule, is(i));
		}
	}
	
	@Test
	public void onePosition() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		// check the pos
		assertThat(confSpace.positions.size(), is(1));
		SimpleConfSpace.Position pos = confSpace.positions.get(0);
		assertThat(pos.index, is(0));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("A2"));
	}
	
	@Test
	public void twoPositions() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType);
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		assertThat(confSpace.positions.size(), is(2));
		
		SimpleConfSpace.Position pos = confSpace.positions.get(0);
		assertThat(pos.index, is(0));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("A2"));
		
		pos = confSpace.positions.get(1);
		assertThat(pos.index, is(1));
		assertThat(pos.strand, is(strand));
		assertThat(pos.resNum, is("A42"));
	}
	
	@Test
	public void onePositionShell() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A42").setLibraryRotamers(Strand.WildType);
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		assertThat(strand.mol.getResByPDBResNumber("A2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.dofBounds.size(), is(0));
	}
	
	@Test
	public void oneDiscreteGlycine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers("GLY");
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.dofBounds.size(), is(0));
	}
	
	@Test
	public void oneDiscreteWildTypeAndGlycine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "GLY");
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.dofBounds.size(), is(0));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("GLY"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
	}
	
	@Test
	public void oneDiscreteValine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers("VAL");
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0)); // valine has one chi angle, but this is discrete flex
		assertThat(pmol.dofBounds.size(), is(0));
		
		double expChi1 = template.getRotamericDihedrals(0, 0);
		double obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("A2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
		
		expChi1 = template.getRotamericDihedrals(1, 0);
		obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("A2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
		
		expChi1 = template.getRotamericDihedrals(2, 0);
		obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("A2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));
	}
	
	@Test
	public void oneContinuousValine() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers("VAL").setContinuous(3);
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		ResidueTemplate template = strand.templateLib.getTemplate("VAL");
		
		// check DOFs and bounds
		int[] conf = { 0 };
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		Residue res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1)); // valine has one chi angle
		FreeDihedral dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		assertThat(pmol.dofBounds.size(), is(1));
		double chi1 = template.getRotamericDihedrals(0, 0);
		assertThat(pmol.dofBounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(pmol.dofBounds.getMax(0), isRelatively(chi1 + 3));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		assertThat(pmol.dofBounds.size(), is(1));
		chi1 = template.getRotamericDihedrals(1, 0);
		assertThat(pmol.dofBounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(pmol.dofBounds.getMax(0), isRelatively(chi1 + 3));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));
		
		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral)pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));
		
		assertThat(pmol.dofBounds.size(), is(1));
		chi1 = template.getRotamericDihedrals(2, 0);
		assertThat(pmol.dofBounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(pmol.dofBounds.getMax(0), isRelatively(chi1 + 3));
	}
	
	@Test
	public void wildTypeRotamersOnly() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").addWildTypeRotamers();
		assertThat(strand.mol.getResByPDBResNumber("A2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
	}
	
	@Test
	public void wildTypeAndLibraryRotamers() {
		
		Strand strand = makeStrand();
		strand.flexibility.get("A2").setLibraryRotamers("VAL").addWildTypeRotamers();
		assertThat(strand.mol.getResByPDBResNumber("A2").fullName, is("ALA A   2"));
		
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
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
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
		
		conf = new int[] { 1 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
		
		conf = new int[] { 2 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
		
		conf = new int[] { 3 };
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));
	}
	
	@Test
	public void twoStrands() {
		
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A3").build();
		strand1.flexibility.get("A2").setLibraryRotamers("GLY");
		strand1.flexibility.get("A3").setLibraryRotamers("GLY");
		Strand strand2 = new Strand.Builder(mol).setResidues("A10", "A11").build();
		strand2.flexibility.get("A10").setLibraryRotamers("GLY");
		strand2.flexibility.get("A11").setLibraryRotamers("GLY");
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(strand1, strand2).build();
		
		assertThat(confSpace.positions.size(), is(4));
		
		Position pos = confSpace.positions.get(0);
		assertThat(pos.index, is(0));
		assertThat(pos.strand, is(strand1));
		assertThat(pos.resNum, is("A2"));
		assertThat(pos.resConfs.size(), is(1));
		assertThat(pos.resConfs.get(0).template, is(strand1.templateLib.getTemplate("GLY")));
		
		pos = confSpace.positions.get(1);
		assertThat(pos.index, is(1));
		assertThat(pos.strand, is(strand1));
		assertThat(pos.resNum, is("A3"));
		assertThat(pos.resConfs.size(), is(1));
		assertThat(pos.resConfs.get(0).template, is(strand1.templateLib.getTemplate("GLY")));
		
		pos = confSpace.positions.get(2);
		assertThat(pos.index, is(2));
		assertThat(pos.strand, is(strand2));
		assertThat(pos.resNum, is("A10"));
		assertThat(pos.resConfs.size(), is(1));
		assertThat(pos.resConfs.get(0).template, is(strand2.templateLib.getTemplate("GLY")));
		
		pos = confSpace.positions.get(3);
		assertThat(pos.index, is(3));
		assertThat(pos.strand, is(strand2));
		assertThat(pos.resNum, is("A11"));
		assertThat(pos.resConfs.size(), is(1));
		assertThat(pos.resConfs.get(0).template, is(strand2.templateLib.getTemplate("GLY")));
	}
	
	@Test
	public void twoStrandsMakeMolecule() {
		
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A3").build();
		strand1.flexibility.get("A2").setLibraryRotamers("GLY");
		strand1.flexibility.get("A3").setLibraryRotamers("GLY");
		Strand strand2 = new Strand.Builder(mol).setResidues("A10", "A11").build();
		strand2.flexibility.get("A10").setLibraryRotamers("GLY");
		strand2.flexibility.get("A11").setLibraryRotamers("GLY");
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder().addStrands(strand1, strand2).build();
		
		ParametricMolecule pmol = confSpace.makeMolecule(new int[] { 0, 0, 0, 0 });
		
		// all the residues should be there
		assertThat(pmol.mol.residues.size(), is(4));
		assertThat(pmol.mol.residues.get(0).getPDBResNumber(), is("A2"));
		assertThat(pmol.mol.residues.get(1).getPDBResNumber(), is("A3"));
		assertThat(pmol.mol.residues.get(2).getPDBResNumber(), is("A10"));
		assertThat(pmol.mol.residues.get(3).getPDBResNumber(), is("A11"));
		
		// and the indices should be set
		for (int i=0; i<pmol.mol.residues.size(); i++) {
			assertThat(pmol.mol.residues.get(i).indexInMolecule, is(i));
		}
	}

	@Test
	public void twoStrandsShell() {

		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers("GLY");
		strand1.flexibility.get("A3").setLibraryRotamers("GLY");
		Strand strand2 = new Strand.Builder(mol).setResidues("A11", "A20").build();
		strand2.flexibility.get("A11").setLibraryRotamers("GLY");
		strand2.flexibility.get("A12").setLibraryRotamers("GLY");
		SimpleConfSpace separateConfSpace = new SimpleConfSpace.Builder()
			.addStrands(strand1, strand2)
			.setShellDistance(9)
			.build();

		Strand combinedStrand = new Strand.Builder(mol).setResidues("A2", "A20").build();
		combinedStrand.flexibility.get("A2").setLibraryRotamers("GLY");
		combinedStrand.flexibility.get("A3").setLibraryRotamers("GLY");
		combinedStrand.flexibility.get("A11").setLibraryRotamers("GLY");
		combinedStrand.flexibility.get("A12").setLibraryRotamers("GLY");
		SimpleConfSpace combinedConfSpace = new SimpleConfSpace.Builder()
			.addStrands(combinedStrand)
			.setShellDistance(9)
			.build();

		assertThat(separateConfSpace.shellResNumbers, is(combinedConfSpace.shellResNumbers));
	}
}
