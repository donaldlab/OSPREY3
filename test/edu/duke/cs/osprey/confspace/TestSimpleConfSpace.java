/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.confspace;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunctionFactory;
import edu.duke.cs.osprey.markstar.MARKStar;
import edu.duke.cs.osprey.markstar.framework.MARKStarBound;
import edu.duke.cs.osprey.parallelism.Parallelism;
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
		int[] conf = {0};
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
		for (int i = 0; i < pmol.mol.residues.size(); i++) {
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
		int[] conf = {0};
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
		int[] conf = {0};
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
		int[] conf = {0};
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("2").template.name, is("ALA"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));

		conf = new int[]{1};
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
		int[] conf = {0};
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0)); // valine has one chi angle, but this is discrete flex
		assertThat(pmol.dofBounds.size(), is(0));

		double expChi1 = template.getRotamericDihedrals(0, 0);
		double obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("A2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));

		conf = new int[]{1};
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));

		expChi1 = template.getRotamericDihedrals(1, 0);
		obsChi1 = new FreeDihedral(pmol.mol.getResByPDBResNumber("A2"), 0).measureDihedralDegrees();
		assertThat(obsChi1, isRelatively(expChi1));

		conf = new int[]{2};
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
		int[] conf = {0};
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		Residue res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));

		assertThat(pmol.dofs.size(), is(1)); // valine has one chi angle
		FreeDihedral dof = (FreeDihedral) pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));

		assertThat(pmol.dofBounds.size(), is(1));
		double chi1 = template.getRotamericDihedrals(0, 0);
		assertThat(pmol.dofBounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(pmol.dofBounds.getMax(0), isRelatively(chi1 + 3));

		conf = new int[]{1};
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));

		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral) pmol.dofs.get(0);
		assertThat(dof.getResidue(), is(res));
		assertThat(dof.getDihedralNumber(), is(0));

		assertThat(pmol.dofBounds.size(), is(1));
		chi1 = template.getRotamericDihedrals(1, 0);
		assertThat(pmol.dofBounds.getMin(0), isRelatively(chi1 - 3));
		assertThat(pmol.dofBounds.getMax(0), isRelatively(chi1 + 3));

		conf = new int[]{2};
		pmol = confSpace.makeMolecule(conf);
		res = pmol.mol.getResByPDBResNumber("A2");
		assertThat(res.template.name, is("VAL"));

		assertThat(pmol.dofs.size(), is(1));
		dof = (FreeDihedral) pmol.dofs.get(0);
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
		int[] conf = {0};
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
		int[] conf = {0};
		ParametricMolecule pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));

		conf = new int[]{1};
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));

		conf = new int[]{2};
		pmol = confSpace.makeMolecule(conf);
		assertThat(pmol.mol.getResByPDBResNumber("A2").template.name, is("VAL"));
		assertThat(pmol.dofs.size(), is(0));
		assertThat(pmol.dofBounds.size(), is(0));

		conf = new int[]{3};
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

		ParametricMolecule pmol = confSpace.makeMolecule(new int[]{0, 0, 0, 0});

		// all the residues should be there
		assertThat(pmol.mol.residues.size(), is(4));
		assertThat(pmol.mol.residues.get(0).getPDBResNumber(), is("A2"));
		assertThat(pmol.mol.residues.get(1).getPDBResNumber(), is("A3"));
		assertThat(pmol.mol.residues.get(2).getPDBResNumber(), is("A10"));
		assertThat(pmol.mol.residues.get(3).getPDBResNumber(), is("A11"));

		// and the indices should be set
		for (int i = 0; i < pmol.mol.residues.size(); i++) {
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

	/**
	 * Make a copy of a mutable confspace that excludes flexible residues.
	 *
	 * Test to make sure that these residues are indeed gone from the confspace
	 */
	@Test
	public void testFlexibleCopy() {
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "LYS");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "VAL");

		SimpleConfSpace mutableConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(9)
				.build();

		SimpleConfSpace flexibleConfSpace = mutableConfSpace.makeFlexibleCopy();

		assertThat(flexibleConfSpace.immutablePositions, is(mutableConfSpace.immutablePositions));
		assertThat(flexibleConfSpace.shellResNumbers, is(mutableConfSpace.shellResNumbers));
		assertThat(flexibleConfSpace.mutablePositions.isEmpty(), is(true));
	}

	/**
	 * Make a copy of a mutable confspace that excludes specific residues.
	 * In this case, we exclude one mutable, one flexible, and one static.
	 *
	 * Test to make sure that these residues are indeed gone from the confspace
	 */
	@Test
	public void testExcludedResNumCopy() {
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "LYS");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "VAL");

		SimpleConfSpace originalConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		Set<String> toExclude = new HashSet<>();
		toExclude.add("A2");
		toExclude.add("A4");
		toExclude.add("A10");

		SimpleConfSpace excludedConfSpace = originalConfSpace.makeCopyExcludingResNums(toExclude);

		// Assert that the shellResNumbers are properly excluded
		assertThat(originalConfSpace.shellResNumbers.contains("A10"), is(true));
		assertThat(excludedConfSpace.shellResNumbers.contains("A10"), is(false));
		assertThat(originalConfSpace.shellResNumbers, contains("A10", "A7", "A8", "A9"));
		assertThat(excludedConfSpace.shellResNumbers, contains("A7", "A8", "A9"));

		// Assert that the positions are properly excluded
		assertThat(originalConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A2"),
				is(true));
		assertThat(originalConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A4"),
				is(true));
		assertThat(originalConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()),
				contains("A2", "A3", "A4", "A5", "A6"));
		assertThat(excludedConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A2"),
				is(false));
		assertThat(excludedConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A4"),
				is(false));
		assertThat(excludedConfSpace.positions.stream().map(e -> (e.resNum)).collect(Collectors.toList()),
				contains("A3", "A5", "A6"));

		// Assert things about the mutable positions
		assertThat(originalConfSpace.mutablePositions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A2"),
				is(true));
		assertThat(excludedConfSpace.mutablePositions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A2"),
				is(false));

		// Assert things about the flexible positions
		assertThat(originalConfSpace.immutablePositions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A4"),
				is(true));
		assertThat(excludedConfSpace.immutablePositions.stream().map(e -> (e.resNum)).collect(Collectors.toList()).contains("A4"),
				is(false));
	}

	/**
	 * Make a flexible copy of a flexible confspace.
	 * These should be the same confspace, so they should have the same positions and static residues
	 *
	 */
	@Test
	public void testFlexibleCopyOfFlexibleConfSpace(){
		// Making a confspace with flexible residues only
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType);

		SimpleConfSpace originalConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		// Making a confspace copy that removes the mutable residues only
		SimpleConfSpace flexibleConfSpace = originalConfSpace.makeFlexibleCopy();

		assertThat(flexibleConfSpace.immutablePositions, is(originalConfSpace.immutablePositions));
		assertThat(flexibleConfSpace.mutablePositions, is(originalConfSpace.mutablePositions));
		assertThat(flexibleConfSpace.shellResNumbers, is(originalConfSpace.shellResNumbers));
		assertThat(flexibleConfSpace.positions, is(originalConfSpace.positions));
		assertThat(flexibleConfSpace.seqSpace, is(originalConfSpace.seqSpace));
	}

	/**
	 * Make a flexible copy of a flexible confspace.
	 * These should be the same confspace, so their partition functions should be the same
	 */
	@Test
	public void testFlexibleCopyOfFlexibleConfSpacePfuncs(){
		// Making a confspace with flexible residues only
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType);

		SimpleConfSpace originalConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		// Making a confspace copy that removes the mutable residues only
		SimpleConfSpace flexibleConfSpace = originalConfSpace.makeFlexibleCopy();

		PartitionFunction flexPfunc = makeMARKStarPfuncForConfSpace(flexibleConfSpace, 0.68);
		PartitionFunction origPfunc = makeMARKStarPfuncForConfSpace(originalConfSpace, 0.68);

		flexPfunc.compute();
		origPfunc.compute();

		double UBdiff = flexPfunc.getValues().calcUpperBound()
				.divide(origPfunc.getValues().calcUpperBound(), RoundingMode.HALF_UP)
				.subtract(BigDecimal.valueOf(1.0))
				.abs().doubleValue();
		double LBdiff = flexPfunc.getValues().calcLowerBound()
				.divide(origPfunc.getValues().calcLowerBound(), RoundingMode.HALF_UP)
				.subtract(BigDecimal.valueOf(1.0))
				.abs().doubleValue();

		assertThat(flexPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(origPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(UBdiff, lessThan(1e-10));
		assertThat(LBdiff, lessThan(1e-10));
	}

	/**
	 * Make a flexible copy of a mutable confspace, and make an excluded copy excluding only mutable residues
	 * These confspaces should be the same, so they should have the same positions and static residues
	 */
	@Test
	public void testFlexibleCopyAgainstManualCopy() {
		// Making a confspace with both mutable and flexible residues
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "LYS");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "VAL");

		SimpleConfSpace mutableConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		// Making a confspace copy that removes the mutable residues only
		SimpleConfSpace flexibleConfSpace = mutableConfSpace.makeFlexibleCopy();

		// Manually making a confspace that should be the same as the flexible copy, if I'm right
		Set<String> toExclude = new HashSet<>();
		toExclude.add("A2");
		toExclude.add("A3");
		toExclude.add("A6");

		SimpleConfSpace manualFlexibleConfSpace = mutableConfSpace.makeCopyExcludingResNums(toExclude);

		assertThat(flexibleConfSpace.immutablePositions, is(manualFlexibleConfSpace.immutablePositions));
		assertThat(flexibleConfSpace.mutablePositions, is(manualFlexibleConfSpace.mutablePositions));
		assertThat(flexibleConfSpace.shellResNumbers, is(manualFlexibleConfSpace.shellResNumbers));
		assertThat(flexibleConfSpace.positions, is(manualFlexibleConfSpace.positions));
		assertThat(flexibleConfSpace.seqSpace, is(manualFlexibleConfSpace.seqSpace));
	}

	/**
	 * Make a flexible copy of a mutable confspace.
	 *
	 * Ensure that we can calculate the partition function as we would expect.
	 * This suggests that the confspace is valid
	 */
	@Test
	public void testFlexibleCopyPfunc(){
		// Making a confspace with both mutable and flexible residues
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "LYS");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "VAL");

		SimpleConfSpace mutableConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		// Making a confspace copy that removes the mutable residues only
		SimpleConfSpace flexibleConfSpace = mutableConfSpace.makeFlexibleCopy();

		PartitionFunction flexPfunc = makeMARKStarPfuncForConfSpace(flexibleConfSpace, 0.68);
		flexPfunc.compute();

		assertThat(flexPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
	}

	/**
	 * Make a flexible copy of a mutable confspace, and make an excluded copy excluding only mutable residues
	 * These confspaces should be the same, so they should return the same pfunc value
	 */
	@Test
	public void testFlexibleCopyVsExcludedCopyPfunc() {
		// Making a confspace with both mutable and flexible residues
		Strand strand1 = new Strand.Builder(mol).setResidues("A2", "A10").build();
		strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG");
		strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "LYS");
		strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType);
		strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType, "VAL");

		SimpleConfSpace mutableConfSpace = new SimpleConfSpace.Builder()
				.addStrands(strand1)
				.setShellDistance(100)
				.build();

		// Making a confspace copy that removes the mutable residues only
		SimpleConfSpace flexibleConfSpace = mutableConfSpace.makeFlexibleCopy();

		// Manually making a confspace that should be the same as the flexible copy, if I'm right
		Set<String> toExclude = new HashSet<>();
		toExclude.add("A2");
		toExclude.add("A3");
		toExclude.add("A6");

		SimpleConfSpace manualFlexibleConfSpace = mutableConfSpace.makeCopyExcludingResNums(toExclude);

		PartitionFunction flexPfunc = makeMARKStarPfuncForConfSpace(flexibleConfSpace, 0.68);
		PartitionFunction manualPfunc = makeMARKStarPfuncForConfSpace(manualFlexibleConfSpace, 0.68);

        flexPfunc.compute();
        manualPfunc.compute();

        double UBdiff = flexPfunc.getValues().calcUpperBound()
				.divide(manualPfunc.getValues().calcUpperBound(), RoundingMode.HALF_UP)
				.subtract(BigDecimal.valueOf(1.0))
				.abs().doubleValue();
		double LBdiff = flexPfunc.getValues().calcLowerBound()
				.divide(manualPfunc.getValues().calcLowerBound(), RoundingMode.HALF_UP)
				.subtract(BigDecimal.valueOf(1.0))
				.abs().doubleValue();

        assertThat(flexPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(manualPfunc.getStatus(), is(PartitionFunction.Status.Estimated));
		assertThat(UBdiff, lessThan(1e-10));
		assertThat(LBdiff, lessThan(1e-10));
	}

	/**
	 * Computes a single partition function for the wild-type sequence only
	 */
	private PartitionFunction makeMARKStarPfuncForConfSpace(SimpleConfSpace confSpace, double epsilon){
		// Set up partition function requirements
		Parallelism parallelism = Parallelism.makeCpu(4);
		ForcefieldParams ffparams = new ForcefieldParams();

		// how should we compute energies of molecules?
		try (EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
				.setParallelism(parallelism)
				.build()) {
			// how should we define energies of conformations?
			ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
							.build()
							.calcReferenceEnergies()
					)
					.build();

			// BBK* needs rigid energies too
			EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
					.setIsMinimizing(false)
					.build();
			ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);
			PartitionFunctionFactory pfuncFactory = new PartitionFunctionFactory(confSpace, confEcalcMinimized, "pfunc");
			pfuncFactory.setUseMARKStar(confEcalcRigid);
			// filter the global sequence to this conf space
			Sequence sequence = confSpace.seqSpace.makeWildTypeSequence();
			// make the partition function
			RCs rcs = sequence.makeRCs(confSpace);

			PartitionFunction pfunc = pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);
			return pfunc;
		}
	}
}
