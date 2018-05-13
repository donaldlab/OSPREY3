/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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
