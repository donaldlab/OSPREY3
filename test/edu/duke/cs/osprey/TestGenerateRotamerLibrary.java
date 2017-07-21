package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.Molecule;
import org.junit.Test;

import edu.duke.cs.osprey.restypes.ResidueTemplate;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

public class TestGenerateRotamerLibrary {

	@Test
	public void test()
	throws Exception {

		Molecule mol = PDBIO.readFile("examples/4NPD/4NPD.pdb");
		ResidueTemplateLibrary templateLib = new ResidueTemplateLibrary.Builder(ForcefieldParams.Forcefield.AMBER)
			.addMoleculeForWildTypeRotamers(mol)
			.build();

		// match a strand to the lib
		Strand strand = new Strand.Builder(mol)
			.setTemplateLibrary(templateLib)
			.build();

		// check the templates and rotamers
		for (Residue res : strand.mol.residues) {

			ResidueTemplate templ = templateLib.getOrMakeWildTypeTemplate(res);

			assertThat(templ.name, is(res.template.name));
			assertThat(templ.getNumRotamers(0, 0), is(1 + strand.mol.getAlternates(res.indexInMolecule).size()));
		}
	}
}
