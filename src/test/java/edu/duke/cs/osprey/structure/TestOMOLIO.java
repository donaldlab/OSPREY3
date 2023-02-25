package edu.duke.cs.osprey.structure;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;

import edu.duke.cs.osprey.confspace.Strand;
import org.junit.jupiter.api.Test;

import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;


public class TestOMOLIO {

	@Test
	public void test() {

		Molecule mol1 = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build().mol;
		mol1.name = "mol";

		String toml = OMOLIO.write(mol1);

		Molecule mol2 = OMOLIO.read(toml);

		// make sure mol1 and mol2 are the same
		assertThat(mol2.name, is(mol1.name));
		assertThat(mol2.residues.size(), is(mol1.residues.size()));
		for (int r=0; r<mol1.residues.size(); r++) {
			Residue res1 = mol1.residues.get(r);
			Residue res2 = mol2.residues.get(r);

			assertThat(res2.fullName, is(res1.fullName));

			assertThat(res2.atoms.size(), is(res1.atoms.size()));
			for (Atom atom1 : res1.atoms) {
				Atom atom2 = res2.getAtomByName(atom1.name);
				assertThat(atom2, is(not(nullValue())));

				assertThat(atom2.elementType, is(atom1.elementType));
				assertThat(atom2.getCoords(), is(atom1.getCoords()));
				assertThat(atom2.bonds.size(), is(atom1.bonds.size()));

				// check the bonds
				Function<Atom, Set<String>> dumpBonds = atom -> atom.bonds.stream()
					.map(bonded -> bonded.res.getPDBResNumber() + "_" + bonded.name)
					.collect(Collectors.toSet());
				assertThat(dumpBonds.apply(atom2), is(dumpBonds.apply(atom1)));
			}
		}
	}
}
