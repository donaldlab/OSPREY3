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
